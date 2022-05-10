import matplotlib.pyplot as plt
import matplotlib as mpl


# style
mpl.rcParams['text.usetex'] = True
mpl.rcParams['text.latex.preamble'] = \
        r"\usepackage[scaled]{helvet} \renewcommand{\familydefault}{\sfdefault} \usepackage[T1]{fontenc} \renewcommand{\rmfamily}{\scriptsize}"
mpl.rcParams['mathtext.fontset'] = 'custom'
mpl.rcParams['mathtext.it'] = 'Helvetica:italic'
plt.rc('font', family='Helvetica', size=6)
plt.rc('lines', linewidth=0.5, markersize=1.5)
plt.rc('axes', linewidth=0.5, xmargin=0, ymargin=0, labelpad=1)

plt.rc('xtick', direction='inout')
plt.rc('xtick.major', size=2.5, width=0.5)
plt.rc('xtick.minor', size=1, width=0.5)

plt.rc('ytick', direction='inout')
plt.rc('ytick.major', size=2.5, width=0.5)
plt.rc('ytick.minor', size=1, width=0.5)

plt.rcParams['lines.dashed_pattern'] = [2, 4]
plt.rcParams['figure.figsize'] = (1.0, 1.0)  # in inches
plt.rcParams['figure.dpi'] = 1200
plt.rcParams["axes.formatter.use_mathtext"] = True  # 10^x instead of 1ex
plt.rcParams["figure.facecolor"] = (1.0, 1.0, 1.0, 0)  # white with alpha = 0%
plt.rcParams["savefig.facecolor"] = (1.0, 1.0, 1.0, 0)  # white with alpha = 0%
# end style

# common colors

# -- infected receiver
SND_BLUE =       '#03b5f0ff'
SND_BLUE_LIGHT = '#a0e5ffff'

# -- receiver
# ECOLI_PINK =       '#ff8080ff'  # border of E.coli image
# ECOLI_PINK_LIGHT = '#FFB3B3ff'  # middle of E.coli images: '#ffd5d5ff'
ECOLI_PINK =       '#cd5b45ff'
ECOLI_PINK_LIGHT = '#FFB3B3ff'  # middle of E.coli images: '#ffd5d5ff'

# -- infected receiver
RCV_BLUE =       '#248dae'
RCV_BLUE_LIGHT = '#53badb'
RCV_INFECTED   = '#e41ce3ff' # same as plasmid color

# -- phage
# PHAGE_GREEN =       '#87c57aff'
# PHAGE_GREEN_LIGHT = '#c6e6cdff'
PHAGE_GREEN =       '#006400ff'
PHAGE_GREEN_LIGHT = '#c6e6cdff'

# -- snd + rcv
BOTH_ORANGE = '#CC8400'
BOTH_ORANGE_LIGHT = '#FFD17D'

# nice other colors
GREEN_DARK =  "#1b7019"
GREEN_LIGHT = "#50c24c"

RED_DARK =  "#b35744"
RED_LIGHT = "#fc9c88"

GRAY_LIGHT = "#9c9c9c"

BLUE_DARK  = '#1f77b4'
BLUE_LIGHT = '#53badb'