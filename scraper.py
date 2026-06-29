import re
import uuid
from pathlib import Path
from typing import List, Dict, Any, Optional, Tuple

try:
    import fitz
except ImportError:
    raise SystemExit("Установите PyMuPDF: pip install PyMuPDF")

import pandas as pd

# ------------------------------------------------------------
# Пути и настройки
# ------------------------------------------------------------
ROOT = Path(__file__).resolve().parent
PDF_FOLDER = ROOT / "pdfs"
OUT_CSV = ROOT / "chem_data.csv"

# Минимальный набор колонок
COLUMNS = [
    "record_id", "source_pdf", "source_doi", "source_page", "source_table",
    "compound_id", "endpoint_raw", "value_raw", "unit_raw", "qualifier",
    "target_raw", "target_std"
]

TARGET_DICT = {
    "cb1": "CNR1", "cb2": "CNR2",
    "mor": "OPRM1", "mop": "OPRM1", "mu": "OPRM1", "μ": "OPRM1",
    "kappa": "OPRK1", "κ": "OPRK1", "kor": "OPRK1",
    "delta": "OPRD1", "δ": "OPRD1", "dor": "OPRD1",
    "nop": "OPRL1", "orl1": "OPRL1",
    "npsr": "NPSR",
    "herg": "KCNH2",
    "l": "OPRM1", "d": "OPRD1", "j": "OPRK1",
}


# ------------------------------------------------------------
# Базовые утилиты
# ------------------------------------------------------------

def new_record(**fields) -> Dict[str, Any]:
    rec = {c: "" for c in COLUMNS}
    rec["record_id"] = fields.pop("record_id", str(uuid.uuid4()))
    rec.update(fields)
    return rec


def extract_meta(pdf_path: Path) -> Dict[str, str]:
    doc = fitz.open(pdf_path)
    meta = doc.metadata or {}
    doc.close()
    doi = ""
    title = meta.get("title", "")
    for candidate in (title, meta.get("subject", "")):
        m = re.search(r"10\.\d{4,9}/[-._;()/:A-Z0-9]+", candidate or "", re.I)
        if m:
            doi = m.group(0)
            break
    return {"doi": doi, "title": title}


def read_all_text(pdf_path: Path) -> str:
    doc = fitz.open(pdf_path)
    text = "\n".join(page.get_text() for page in doc)
    doc.close()
    return text


# ------------------------------------------------------------
# Извлечение таблиц: camelot + pdfplumber
# ------------------------------------------------------------
def fetch_tables(pdf_path: Path) -> List[Dict]:
    out = []
    try:
        import camelot
        tables = camelot.read_pdf(str(pdf_path), pages="all", flavor="lattice")
        for table in tables:
            rows = table.df.values.tolist()
            out.append({"page": table.page, "rows": rows})
    except Exception:
        pass

    if not out:
        try:
            import pdfplumber
            with pdfplumber.open(pdf_path) as pdf:
                for i, page in enumerate(pdf.pages, 1):
                    for table in page.extract_tables() or []:
                        if table and len(table) > 1:
                            out.append({"page": i, "rows": table})
        except ImportError:
            pass
    return out


def clean_cell(cell) -> str:
    if cell is None:
        return ""
    return re.sub(r"\s+", " ", str(cell).replace("\n", " ")).strip()


def parse_measurement(raw: str) -> Tuple[Optional[float], Optional[str], Optional[str]]:
    if not raw or not raw.strip():
        return (None, None, None)
    s = raw.strip()
    op_match = re.match(r'^([><=~])\s*', s)
    operator = op_match.group(1) if op_match else '='
    if op_match:
        s = s[op_match.end():].strip()
    s = re.sub(r'\([^)]*\)', '', s)
    s = re.sub(r'±\s*[\d.]+', '', s)
    num_match = re.search(r'^([\d.]+(?:[eE][+-]?\d+)?)', s)
    if not num_match:
        return (None, None, operator)
    num = float(num_match.group(1))
    unit_part = s[num_match.end():].strip()
    unit_match = re.search(r'^([µμu]?M|nM|μM|uM|mM|pM|fM|M|pX|log)', unit_part, re.I)
    if unit_match:
        unit = unit_match.group(1)
        unit = unit.upper().replace("µ", "u").replace("μ", "u")
        if unit in ["UM", "U M"]:
            unit = "μM"
        elif unit in ["NM", "N M"]:
            unit = "nM"
        elif unit in ["MM", "M M"]:
            unit = "mM"
        elif unit in ["PM", "P M"]:
            unit = "pM"
        elif unit in ["FM", "F M"]:
            unit = "fM"
        elif unit in ["PX", "P X"]:
            unit = "pX"
        else:
            unit = unit
    else:
        unit = unit_part or None
    return (num, unit, operator)


def standardize_target(raw: str) -> str:
    if not raw:
        return ""
    key = raw.lower().strip()
    if key in TARGET_DICT:
        return TARGET_DICT[key]
    for suffix in (" receptor", " binding", "binding", "receptor"):
        if key.endswith(suffix):
            sub = key[:-len(suffix)].strip()
            if sub in TARGET_DICT:
                return TARGET_DICT[sub]
    return key.upper() if key.isascii() and len(key) <= 6 else key


# ------------------------------------------------------------
# Парсеры
# ------------------------------------------------------------
def parse_fulton2008(pdf_path: Path) -> List[Dict]:
    records = []
    meta = extract_meta(pdf_path)
    for tbl in fetch_tables(pdf_path):
        merged = " ".join(clean_cell(c) for row in tbl["rows"] for c in row)
        if "Ki" not in merged:
            continue
        for m in re.finditer(
                r"(\d+[a-z](?:\([^)]+\))?|[A-Z]+-\d+)\s+"
                r"([\d.]+)(?:[±+\-]\s*([\d.]+))?\s+"
                r"([\d.]+)(?:[±+\-]\s*([\d.]+))?\s+"
                r"([\d.]+)(?:[±+\-]\s*([\d.]+))?",
                merged,
        ):
            comp, ki_l, sem_l, ki_d, sem_d, ki_k, sem_k = m.groups()
            for target, val in (("l", ki_l), ("d", ki_d), ("j", ki_k)):
                if not val:
                    continue
                num, unit, op = parse_measurement(val)
                if num is None:
                    continue
                records.append(new_record(
                    source_pdf=pdf_path.name, source_doi=meta["doi"],
                    source_page=tbl["page"], source_table="Table 1",
                    compound_id=comp,
                    endpoint_raw="Ki",
                    value_raw=num, unit_raw=unit, qualifier=op,
                    target_raw=target, target_std=standardize_target(target)
                ))
    # fallback
    if not records:
        text = read_all_text(pdf_path)
        # ищем любые числа с nM
        for m in re.finditer(r"(MCL-\d+)\s+.*?([\d.]+)\s*nM", text):
            comp, val = m.groups()
            num, unit, op = parse_measurement(val)
            if num is None:
                continue
            records.append(new_record(
                source_pdf=pdf_path.name, source_doi=meta["doi"],
                source_page="", source_table="prose",
                compound_id=comp,
                endpoint_raw="Ki",
                value_raw=num, unit_raw="nM", qualifier=op,
                target_raw="l", target_std="OPRM1"
            ))
    return records


def parse_yang2009(pdf_path: Path) -> List[Dict]:
    records = []
    meta = extract_meta(pdf_path)
    for tbl in fetch_tables(pdf_path):
        merged = " ".join(clean_cell(c) for row in tbl["rows"] for c in row)
        compact = merged.replace(" ", "")
        if "NOPKi" not in compact and "NOP Ki" not in merged:
            continue
        ranges = [(14, 25, 1), (41, 47, 2), (26, 32, 3), (33, 40, 4), (47, 50, 5)]
        for min_id, max_id, table_num in ranges:
            for m in re.finditer(
                    r"\b(\d{1,2})\s+([A-Za-z0-9,\-]+(?:\s+[A-Za-z0-9,\-]+)?)\s+([\d.]+|>1000|nd)\s+([\d.]+|>1000|nd)?",
                    merged,
            ):
                comp_id, subst, nop_val, mop_val = m.groups()
                cid = int(comp_id)
                if cid < min_id or cid > max_id:
                    continue
                label = f"{comp_id} ({subst.strip()})"
                for target, val in (("NOP", nop_val), ("MOP", mop_val)):
                    if not val or val == "nd":
                        continue
                    num, unit, op = parse_measurement(val)
                    if num is None:
                        continue
                    records.append(new_record(
                        source_pdf=pdf_path.name, source_doi=meta["doi"],
                        source_page=tbl["page"], source_table=f"Table {table_num}",
                        compound_id=label,
                        endpoint_raw="Ki",
                        value_raw=num, unit_raw=unit, qualifier=op,
                        target_raw=target, target_std=standardize_target(target)
                    ))
    return records


def parse_dolle2009(pdf_path: Path) -> List[Dict]:
    records = []
    meta = extract_meta(pdf_path)
    for tbl in fetch_tables(pdf_path):
        merged = " ".join(clean_cell(c) for row in tbl["rows"] for c in row)
        if "Ki(nM)" not in merged:
            continue
        for m in re.finditer(
                r"\b(\d{1,2})\s+([\d.]+|c)\s+([\d.]+|c)\s+([\d.]+|c)\s+([\d.]+|c)?\s+([\d.]+|c)?\s+([\d.]+|c)?\s+([\d.]+|c)?",
                merged,
        ):
            comp, ki_j, ki_l, ki_d, ic50_j, ic50_l, ic50_d, _ = m.groups()
            if int(comp) > 22:
                continue
            for target, val in (("j", ki_j), ("l", ki_l), ("d", ki_d)):
                if not val or val == "c":
                    continue
                num, unit, op = parse_measurement(val)
                if num is None:
                    continue
                records.append(new_record(
                    source_pdf=pdf_path.name, source_doi=meta["doi"],
                    source_page=tbl["page"], source_table="Table 1",
                    compound_id=comp,
                    endpoint_raw="Ki",
                    value_raw=num, unit_raw=unit, qualifier=op,
                    target_raw=target, target_std=standardize_target(target)
                ))
    if not records:
        text = read_all_text(pdf_path)
        for m in re.finditer(r"compound\s+(\d+)\s+.*?([\d.]+)\s*nM", text):
            comp, val = m.groups()
            num, unit, op = parse_measurement(val)
            if num is None:
                continue
            records.append(new_record(
                source_pdf=pdf_path.name, source_doi=meta["doi"],
                source_page="", source_table="prose",
                compound_id=comp,
                endpoint_raw="Ki",
                value_raw=num, unit_raw="nM", qualifier=op,
                target_raw="l", target_std="OPRM1"
            ))
    return records


def parse_kobayashi2009(pdf_path: Path) -> List[Dict]:
    records = []
    meta = extract_meta(pdf_path)
    for tbl in fetch_tables(pdf_path):
        merged = " ".join(clean_cell(c) for row in tbl["rows"] for c in row)
        if "ORL1" not in merged:
            continue
        for m in re.finditer(
                r"\b(\d{1,2})\s+(?:[^\d]{0,20}\s+)?([\d.]+|>1000)\s+([\d.]+|>1000)?\s*([\d.]+|>1000)?",
                merged,
        ):
            comp, v1, v2, v3 = m.groups()
            if int(comp) > 31:
                continue
            val = v1
            if not val:
                continue
            q = ">" if str(val).startswith(">") else ""
            val_clean = str(val).replace(">1000", "1000").lstrip(">")
            num, unit, op = parse_measurement(val_clean)
            if num is None:
                continue
            records.append(new_record(
                source_pdf=pdf_path.name, source_doi=meta["doi"],
                source_page=tbl["page"], source_table=f"Table (p{tbl['page']})",
                compound_id=comp,
                endpoint_raw="IC50",
                value_raw=num, unit_raw=unit, qualifier=q or op,
                target_raw="ORL1", target_std=standardize_target("ORL1")
            ))
    return records


def parse_iyer2012(pdf_path: Path) -> List[Dict]:
    records = []
    meta = extract_meta(pdf_path)
    text = read_all_text(pdf_path)
    patterns = [
        r"compound\s+(\d+[a-z]?)\s+.*?Ki\s*[=¼]\s*([\d.]+)\s*nM",
        r"(\d+[a-z]?)\s+\(Ki\s*[=¼]\s*([\d.]+)\s*nM\)",
        r"(\d+[a-z]?)\s+.*?Ki\s*[=¼]\s*([\d.]+)\s*nM",
        r"Ki\s*[=¼]\s*([\d.]+)\s*nM.*?(\d+[a-z]?) for the m",
    ]
    seen = set()
    for pat in patterns:
        for m in re.finditer(pat, text, re.I):
            if len(m.groups()) == 2:
                g1, g2 = m.groups()
                if re.match(r"[\d.]+", g1):
                    value, comp = g1, g2
                else:
                    comp, value = g1, g2
            else:
                continue
            key = (comp.lower(), value)
            if key in seen:
                continue
            seen.add(key)
            num, unit, op = parse_measurement(value)
            if num is None:
                continue
            records.append(new_record(
                source_pdf=pdf_path.name, source_doi=meta["doi"],
                source_page="", source_table="prose",
                compound_id=comp,
                endpoint_raw="Ki",
                value_raw=num, unit_raw="nM", qualifier=op,
                target_raw="m-opioid", target_std="OPRM1"
            ))
    return records


def parse_pettersson2009(pdf_path: Path) -> List[Dict]:
    records = []
    meta = extract_meta(pdf_path)
    text = read_all_text(pdf_path)
    # ищем pKi
    for m in re.finditer(r"(12[jk])\s+.*?pKi\s*=\s*([\d.]+)", text):
        comp, val = m.groups()
        num, unit, op = parse_measurement(val)
        if num is None:
            continue
        records.append(new_record(
            source_pdf=pdf_path.name, source_doi=meta["doi"],
            source_page="", source_table="prose",
            compound_id=comp,
            endpoint_raw="pKi",
            value_raw=num, unit_raw="pX", qualifier=op,
            target_raw="CB1", target_std="CNR1"
        ))
    # если ничего не нашли, ищем просто pIC50
    if not records:
        for m in re.finditer(r"(12[a-z])\s+.*?pIC50\s*=\s*([\d.]+)", text):
            comp, val = m.groups()
            num, unit, op = parse_measurement(val)
            if num is None:
                continue
            records.append(new_record(
                source_pdf=pdf_path.name, source_doi=meta["doi"],
                source_page="", source_table="prose",
                compound_id=comp,
                endpoint_raw="pIC50",
                value_raw=num, unit_raw="pX", qualifier=op,
                target_raw="CB1", target_std="CNR1"
            ))
    return records


def parse_guerrini2009(pdf_path: Path) -> List[Dict]:
    records = []
    meta = extract_meta(pdf_path)
    text = read_all_text(pdf_path)
    # ищем pKB
    for m in re.finditer(r"\[tBu-D-Gly5\]NPS.*?pKB\s*=\s*([\d.]+)", text):
        val = m.group(1)
        num, unit, op = parse_measurement(val)
        if num is None:
            continue
        records.append(new_record(
            source_pdf=pdf_path.name, source_doi=meta["doi"],
            source_page="", source_table="prose",
            compound_id="[tBu-D-Gly5]NPS",
            endpoint_raw="pKB",
            value_raw=num, unit_raw="pX", qualifier=op,
            target_raw="NPSR", target_std="NPSR"
        ))
    if not records:
        for m in re.finditer(r"NPS\s+.*?pEC50\s*=\s*([\d.]+)", text):
            val = m.group(1)
            num, unit, op = parse_measurement(val)
            if num is None:
                continue
            records.append(new_record(
                source_pdf=pdf_path.name, source_doi=meta["doi"],
                source_page="", source_table="prose",
                compound_id="NPS",
                endpoint_raw="pEC50",
                value_raw=num, unit_raw="pX", qualifier=op,
                target_raw="NPSR", target_std="NPSR"
            ))
    return records


def parse_naltrexamine(pdf_path: Path) -> List[Dict]:
    records = []
    meta = extract_meta(pdf_path)
    text = read_all_text(pdf_path)
    # ищем Ki в тексте
    for m in re.finditer(r"(\d+[a-z]?)\s+\(Ki\s*=\s*([\d.]+)\s*nM\)", text):
        comp, val = m.groups()
        num, unit, op = parse_measurement(val)
        if num is None:
            continue
        records.append(new_record(
            source_pdf=pdf_path.name, source_doi=meta["doi"],
            source_page="", source_table="prose",
            compound_id=comp,
            endpoint_raw="Ki",
            value_raw=num, unit_raw="nM", qualifier=op,
            target_raw="MOR", target_std="OPRM1"
        ))
    return records


# ------------------------------------------------------------
# Главная функция
# ------------------------------------------------------------
def main():
    if not PDF_FOLDER.exists():
        PDF_FOLDER.mkdir(parents=True, exist_ok=True)
        print(f"Папка {PDF_FOLDER} создана. Поместите PDF-файлы и запустите снова.")
        return

    pdfs = list(PDF_FOLDER.glob("*.pdf"))
    if not pdfs:
        print("В папке pdfs/ нет PDF-файлов.")
        return

    parser_map = {
        "08008214": parse_fulton2008,
        "09003679": parse_yang2009,
        "09006222": parse_dolle2009,
        "09006258": parse_kobayashi2009,
        "02235234": parse_iyer2012,
        "dibenzothiazep": parse_pettersson2009,
        "neuropeptide": parse_guerrini2009,
        "naltrexamine": parse_naltrexamine,
    }

    all_rows = []
    for pdf in pdfs:
        print(f"Обработка {pdf.name}...")
        name = pdf.name.lower()
        parser = None
        for pattern, func in parser_map.items():
            if pattern in name:
                parser = func
                break
        if parser is None:
            print(f"  Нет парсера для {pdf.name}, пропускаем.")
            continue
        records = parser(pdf)
        all_rows.extend(records)
        print(f"  Извлечено {len(records)} записей.")

    if not all_rows:
        print("Не извлечено ни одной записи.")
        return

    df = pd.DataFrame(all_rows)
    for c in COLUMNS:
        if c not in df.columns:
            df[c] = ""
    df = df[COLUMNS]
    df.to_csv(OUT_CSV, index=False, encoding="utf-8-sig")
    print(f"База данных сохранена: {OUT_CSV} ({len(all_rows)} записей)")


if __name__ == "__main__":
    main()
