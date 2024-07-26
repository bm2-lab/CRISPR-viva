import argparse
import sys
from src import ui


prog = 'CRISPR-viva'
usage = 'python crisprviva.py <command> [options]'
p = argparse.ArgumentParser(prog=prog,usage=usage, formatter_class=argparse.RawDescriptionHelpFormatter)
sp = p.add_subparsers(title='Command',metavar='')

p_ref = sp.add_parser('ref', help='Building reference')
ui.parse_ref(p_ref)

p_expr = sp.add_parser('expr', help='RNA Expression profile analysis')
ui.parse_expr(p_expr)

p_elim = sp.add_parser('elim', help='gRNA generation for viral elimination')
ui.parse_elim(p_elim)

p_det = sp.add_parser('det', help='gRNA generation for viral detection')
ui.parse_det(p_det)

p_util = sp.add_parser('util', help='Utilities')
ui.parse_util(p_util)

opts = p.parse_args()

if len(sys.argv) == 1:
    p.print_help()
    sys.exit(1)
elif len(sys.argv) == 2:
    if sys.argv[1] == 'ref':
        p_ref.print_help()
    elif sys.argv[1] == 'expr':
        p_expr.print_help()
    elif sys.argv[1] == 'elim':
        p_elim.print_help()
    elif sys.argv[1] == 'det':
        p_det.print_help()
    elif sys.argv[1] == 'util':
        p_util.print_help()
    else:
        p.print_help()
    sys.exit(1)

if sys.argv[1] == 'ref':
    ui.run_ref(opts)
elif sys.argv[1] == 'expr':
    ui.run_expr(opts)
elif sys.argv[1] == 'elim':
    ui.run_elim(opts)
elif sys.argv[1] == 'det':
    ui.run_det(opts)
elif sys.argv[1] == 'util':
    ui.run_util(opts)