import logging
from lib import log
module_logger = log.setup_logger('root')
module_logger.setLevel(logging.ERROR)


def generate_param_file(set_to, fname_in, fname_out):
    """
    Example:

    set_to = {
        'init_ecoli': ['10', 'float'],
        }
    """
    module_logger.info(f"Generating param file")
    new = ""
    with open(fname_in, 'r') as f:
        for line in f:
            found = False
            for prop in set_to.keys():
                if prop in line.split('\t'):
                    # new line
                    #print(set_to)
                    #print(prop)
                    new += f"{prop}\t\t{set_to[prop][0]}\t\t{set_to[prop][1]}\n"
                    found = True
                    break
            # at the end check if found
            if not found:
                # not changed
                new += line

    with open(fname_out, 'w') as output:
        output.write(new)
    module_logger.info(f"written: {fname_out}")