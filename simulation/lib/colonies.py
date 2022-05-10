def brs_od2platereader_od(od: float):
    """
    convert OD (brs spectro) -> OD (plate reader)

    obtained from script:
    fittingsOD_to_OD/fit.py
    """
    return od * 0.59508123

def platereader_od2brs_od(od: float):
    """
    convert OD (plate reader) -> OD (brs spectro)

    obtained from script:
    fittingsOD_to_OD/fit.py
    """
    return od * 1/0.59508123


def cfu2od(cfu: float, curve='cubic', od_from_brs=True, cell_type='ER2738-Fprime'):

    od2cfu_fun = lambda od: od2cfu(od=od, curve=curve, od_from_brs=od_from_brs, cell_type=cell_type)
    
    # simple binary search
    iterations = 100
    eps = 1
    od_min = 0
    od_max = 10
    for i in range(iterations):
        od_model = (od_min + od_max)/2
        cfu_model = od2cfu_fun(od_model)
        
        if abs(cfu_model - cfu) < eps:
            break

        if cfu_model < cfu:
            od_min = od_model
        else:
            od_max = od_model

    return od_model



def od2cfu(od: float, curve='cubic', od_from_brs=True, cell_type='ER2738-Fprime'):
    """
    obtained from fitting/OD_to_cellcount

    In: OD
    Out: cfu/ml
    """
    if not od_from_brs:
        # convert OD (plate reader) -> OD (brs spectro)
        od = platereader_od2brs_od(od=od)

    if cell_type == 'ER2738-Fprime':
        if od <= 0:
            ret = 0
        elif curve == 'linear':
            ret = 242829272.89818564 * od + 5.645371295905335
        elif curve == 'quad':
            ret = 97408898.69667163 * od**2 + 91073978.31924428 * od + 1.309095379507708
        elif curve == 'cubic':
            ret = 14998665.638420682 * od**3 + 113611796.89482512 * od**2 + -33134.9293802388 * od + -1.6386940602764994
        elif curve == 'exp':
            ret = 18529993.888179254 * exp(1.6596998440406756 * od)
        else:
            assert False, f"unsupported curve: {curve}"

    elif cell_type == 'MG1655':
        if od <= 0:
            ret = 0
        elif curve == 'linear':
            ret = 561423603.3589603 * od + 2.989766383178683
        elif curve == 'quad':
            ret = 231435779.573607 * od**2 + 4.9458022063590725 * od + 698788.3569910543
        elif curve == 'cubic':
            ret = 62752575.54745298 * od**3 + 26124765.864032783 * od**2 + 94474499.37268054 * od + 1.6998478903241843
        elif curve == 'exp':
            ret = 27166797.43821345 * math.exp(1.5223028922167592 * od)
        elif curve == 'factor_cubic':
            # a factor fitted from the ER2738-Fprime cubic
            FACTOR = 1.5
            ret = 14998665.638420682 * od**3 + 113611796.89482512 * od**2 + -33134.9293802388 * od + -1.6386940602764994
            ret *= FACTOR
        else:
            assert False, f"unsupported curve: {curve}"

    else:
        assert False, f"unsupported cell type: {cell_type}"

    return int(ret)


def get_cfu_withminmax(diluted, count,
    fold_per_dilution= 10,
    fraction_plated= 0.1,
    error_model= {'dilution_error': 0.01,
                  'plating_error': 0.01,
                  'count_error': 2}):
    """
    Returns: estimated cfu, [neg error, pos error]
    """
    dilution_error = error_model['dilution_error']
    plating_error = error_model['plating_error']
    count_error = error_model['count_error']
    assert ( count >= count_error )
    measured = count                 * fold_per_dilution                          ** diluted * 1/fraction_plated
    mymin    = (count - count_error) * (fold_per_dilution * (1 - dilution_error)) ** diluted * 1/(fraction_plated * (1+ plating_error))
    mymax    = (count + count_error) * (fold_per_dilution * (1 + dilution_error)) ** diluted * 1/(fraction_plated * (1- plating_error))
    return measured, [[measured - mymin], [mymax - measured]]

