
"""
Скрипт для извлечения данных о биологической активности молекул
из научных PDF-статей. Собирает единую таблицу с данными.
Выход: chem_data.csv
"""

import pdfplumber
import re
import pandas as pd
import os
import glob
from collections import defaultdict


# 1. Вспомогательные функции для извлечения чисел и очистки

def extract_number(cell):
    """Извлекает первое число из строки (поддерживает float и int)."""
    if cell is None:
        return None
    cell_str = str(cell).strip()
    match = re.search(r'([\d.]+)', cell_str)
    if match:
        try:
            return float(match.group(1))
        except ValueError:
            return None
    return None

def clean_text(t):
    """Очищает текст от лишних пробелов и символов."""
    if not t:
        return ''
    return re.sub(r'\s+', ' ', str(t)).strip()


# 2. Парсинг таблиц с активностями (Ki, IC50, EC50)

def parse_activity_table(table, pdf_name):
    """
    Парсит одну таблицу, извлекая названия соединений и значения
    для разных рецепторов. Возвращает список словарей.
    """
    if not table:
        return []

    # Поиск строки-заголовка по ключевым словам
    header_row_idx = None
    for i, row in enumerate(table):
        if row and any(re.search(r'Ki|IC50|EC50|K_i|pKi|pIC50', str(cell), re.IGNORECASE) for cell in row):
            header_row_idx = i
            break
    if header_row_idx is None:
        # Если не нашли, возможно таблица без заголовка – пропускаем
        return []

    headers = [clean_text(cell) for cell in table[header_row_idx]]
    # Определяем индексы столбцов по ключевым словам
    col_indices = {
        'compound': None,
        'mu': None,
        'delta': None,
        'kappa': None,
        'orl1': None,
        'clogp': None,
        'selectivity': None,
        'ec50': None,      # для функциональных тестов
        'emax': None,
        'ke': None
    }

    for idx, h in enumerate(headers):
        h_lower = h.lower()
        if re.search(r'compound|ligand|drug|name|number', h_lower):
            col_indices['compound'] = idx
        if re.search(r'μ|mu|mop|damgo', h_lower):
            col_indices['mu'] = idx
        if re.search(r'δ|delta|dop|naltrindole|dlt', h_lower):
            col_indices['delta'] = idx
        if re.search(r'κ|kappa|kop|u69|u-69', h_lower):
            col_indices['kappa'] = idx
        if re.search(r'orl1|nociceptin|nop', h_lower):
            col_indices['orl1'] = idx
        if re.search(r'clogp|logp|log p', h_lower):
            col_indices['clogp'] = idx
        if re.search(r'selectivity|ratio', h_lower):
            col_indices['selectivity'] = idx
        if re.search(r'ec50|ec 50', h_lower):
            col_indices['ec50'] = idx
        if re.search(r'emax|e max|% stim', h_lower):
            col_indices['emax'] = idx
        if re.search(r'ke|k e', h_lower):
            col_indices['ke'] = idx

    # Если не нашли compound, предположим первый столбец
    if col_indices['compound'] is None:
        col_indices['compound'] = 0

    data = []
    for row in table[header_row_idx+1:]:
        if not row or all(cell is None or clean_text(cell) == '' for cell in row):
            continue

        compound = clean_text(row[col_indices['compound']]) if col_indices['compound'] < len(row) else ''
        if not compound:
            continue

        # Извлекаем числа для каждого рецептора
        mu = extract_number(row[col_indices['mu']]) if col_indices['mu'] is not None and col_indices['mu'] < len(row) else None
        delta = extract_number(row[col_indices['delta']]) if col_indices['delta'] is not None and col_indices['delta'] < len(row) else None
        kappa = extract_number(row[col_indices['kappa']]) if col_indices['kappa'] is not None and col_indices['kappa'] < len(row) else None
        orl1 = extract_number(row[col_indices['orl1']]) if col_indices['orl1'] is not None and col_indices['orl1'] < len(row) else None
        clogp = extract_number(row[col_indices['clogp']]) if col_indices['clogp'] is not None and col_indices['clogp'] < len(row) else None
        ec50 = extract_number(row[col_indices['ec50']]) if col_indices['ec50'] is not None and col_indices['ec50'] < len(row) else None
        emax = extract_number(row[col_indices['emax']]) if col_indices['emax'] is not None and col_indices['emax'] < len(row) else None
        ke = extract_number(row[col_indices['ke']]) if col_indices['ke'] is not None and col_indices['ke'] < len(row) else None

        # Если нет никаких данных по активности – пропускаем
        if all(v is None for v in [mu, delta, kappa, orl1, ec50, ke]):
            continue

        data.append({
            'compound': compound,
            'mu_ki': mu,
            'delta_ki': delta,
            'kappa_ki': kappa,
            'orl1_ki': orl1,
            'clogp': clogp,
            'ec50': ec50,
            'emax': emax,
            'ke': ke,
            'source_pdf': pdf_name
        })
    return data


# 3. Основная функция обработки всех PDF

def process_pdfs(pdf_folder='.'):
    pdf_files = glob.glob(os.path.join(pdf_folder, '*.pdf'))
    if not pdf_files:
        print("PDF-файлы не найдены в папке:", pdf_folder)
        return

    all_activity_data = []

    for pdf_path in pdf_files:
        pdf_name = os.path.basename(pdf_path)
        print(f'Обработка: {pdf_name}')

        # Извлекаем таблицы
        with pdfplumber.open(pdf_path) as pdf:
            for page in pdf.pages:
                tables = page.extract_tables()
                for table in tables:
                    if table and len(table) > 1:
                        parsed = parse_activity_table(table, pdf_name)
                        all_activity_data.extend(parsed)


    # 4. Формирование единой таблицы (денормализованной)

    df_raw = pd.DataFrame(all_activity_data)

    # Удаляем дубликаты по названию и источнику (оставляем первое вхождение)
    df_raw = df_raw.drop_duplicates(subset=['compound', 'source_pdf'], keep='first')

    # Фильтруем: оставляем только записи, где есть данные по хотя бы одному рецептору
    df = df_raw.dropna(subset=['mu_ki', 'delta_ki', 'kappa_ki', 'orl1_ki', 'ec50', 'ke'], how='all')

    # Присваиваем уникальный ID для каждой уникальной молекулы (по названию)
    unique_molecules = df['compound'].unique()
    mol_id_map = {name: f'MOL{str(i+1).zfill(3)}' for i, name in enumerate(unique_molecules)}
    df['molecule_id'] = df['compound'].map(mol_id_map)


    # 5. Создание единого списка записей

    rows = []
    for _, row in df.iterrows():
        mid = row['molecule_id']
        name = row['compound']
        clogp = row['clogp']
        pdf = row['source_pdf']

        # Helper для добавления строки
        def add_activity(target, assay_type, activity_type, value, efficacy=None, selectivity=None):
            if pd.notna(value):
                rows.append({
                    'molecule_id': mid,
                    'name': name,
                    'smiles': '',  # можно заполнить отдельно
                    'clogp': clogp,
                    'target': target,
                    'assay_type': assay_type,
                    'activity_type': activity_type,
                    'value': value,
                    'units': 'nM',
                    'efficacy': efficacy,
                    'selectivity_ratio': selectivity,
                    'source_pdf': pdf
                })

        # Ki для μ
        add_activity('μ opioid receptor', 'binding', 'Ki', row['mu_ki'])
        # Ki для δ
        add_activity('δ opioid receptor', 'binding', 'Ki', row['delta_ki'])
        # Ki для κ
        add_activity('κ opioid receptor', 'binding', 'Ki', row['kappa_ki'])
        # Ki для ORL1
        add_activity('ORL1', 'binding', 'Ki', row['orl1_ki'])
        # EC50 (по умолчанию считаем, что это μ, если не указано иное)
        if pd.notna(row['ec50']):
            add_activity('μ opioid receptor', 'functional', 'EC50', row['ec50'], efficacy=row.get('emax'))
        # Ke (по умолчанию μ)
        if pd.notna(row['ke']):
            add_activity('μ opioid receptor', 'functional', 'Ke', row['ke'])

    df_out = pd.DataFrame(rows)

    # Удаляем дубликаты (могут быть, если одна активность попала дважды)
    df_out = df_out.drop_duplicates()

    # Сортируем для удобства
    df_out = df_out.sort_values(['molecule_id', 'target', 'activity_type'])

    # Сохраняем
    df_out.to_csv('chem_data.csv', index=False, encoding='utf-8')
    print(f"\n✅ Обработка завершена. Создан файл chem_data.csv с {len(df_out)} записями.")

    return df_out


# 6. Запуск

if __name__ == '__main__':
    process_pdfs(pdf_folder='pdfs')
