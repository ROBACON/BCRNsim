import pandas as pd
from statistics import mean

from libs.colonies import od2cfu


def __rename_header(c: str):
    if 'Time' in c:
        if c == 'Time_min':
            return 'Time'
        else:
            return c
    else:
        return c.split(':')[0]


def __get_blank_od(blanks=None, CSV_INPUT_DIR=''):
    if blanks is None:
        return 0

    elif isinstance(blanks, dict):
        fname = blanks['filename']
        df = pd.read_csv(CSV_INPUT_DIR + fname + '.csv')
        return mean([ df[blank].mean() for blank in blanks['fields'] ])

    else:
        return blanks


def odcsv2coloniescsv(FILENAMES,
    CURVE='cubic',
    CSV_INPUT_DIR='',
    cell_type='ER2738-Fprime',
    od_from_brs=True,
    blanks=None):
    """
    example parameters:

    blanks={'filename': 'trace.csv',
            'fields': ['D6','D12','H6','H12']}
    """
    # get blank
    if od_from_brs:
        # blank was automatically set to 0 by the machine
        blank_od = 0
    else:
        blank_od = __get_blank_od(blanks=blanks, CSV_INPUT_DIR=CSV_INPUT_DIR)

    for fname in FILENAMES:
        # read csv
        df = pd.read_csv(CSV_INPUT_DIR + fname + '.csv')
        col_names = list(df.columns)
        #print(col_names)
        col_names_new = [__rename_header(c) for c in col_names]
        #print(col_names_new)
        col_names_new_set = set(col_names_new)
        if not len(col_names_new_set) == len(col_names_new):
            raise AssertionError(f'Generated CSV has not unique headers.\n Headers: {col_names_new}')

        # rename columns
        df.columns = col_names_new

        # calc
        for col_name in col_names_new:
            if 'Time' not in col_name:
                # ATTENTION: assume it is an OD value
                df[col_name] = df.apply( lambda x: od2cfu(od=x[col_name] - blank_od,
                                                          curve=CURVE,
                                                          od_from_brs=od_from_brs,
                                                          cell_type=cell_type),
                                         axis=1 )
 
        # write new file
        out_fname = f"{fname}_OD{CURVE}.csv"
        print(f'--- Generate csv file {out_fname} ---')
        # write comment header
        with open(out_fname, 'w') as out_file:
            out_file.write(f"#--- Do not edit. Generated automatically from {fname} ---\n#{','.join(col_names)}\n")
        # write header and data
        df.to_csv(out_fname, index=False, mode='a')
        print('[done]')
