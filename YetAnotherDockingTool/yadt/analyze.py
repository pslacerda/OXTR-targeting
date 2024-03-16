import re
import sys
from glob import glob
from sklearn.metrics import roc_auc_score, roc_curve
import matplotlib.pyplot as plt
import click


def analyze_vina(filenames, active_re):
    scores = []
    actives = []
    for filename in filenames:
        with open(filename) as file:
            for line in file:
                if line[7:18] == 'VINA RESULT':
                    active = bool(re.search(active_re, filename))
                    score = float(line.split()[3])
                    scores.append(-score)
                    actives.append(active)
                    break
            else:
                raise Exception(f"No score found in '{filename}'.")
    return scores, actives


def analyze_gold(filenames, active_re):
    scores = []
    actives = []
    for filename in filenames:
        with open(filename) as file:
            lines = file.readlines()
            idx = 0
            while idx < len(lines):
                if '@<TRIPOS>MOLECULE' in lines[idx]:
                    idx += 1
                    name = lines[idx].split('|')[0]
                    active = bool(re.search(active_re, name))
                if '> <Gold.Score>' in lines[idx]:
                    idx += 2
                    score = float(lines[idx].split()[0])
                    break
                idx += 1
        actives.append(active)
        scores.append(score)
    return scores, actives


@click.command()
@click.argument('active_regex', nargs=1)
@click.option('-v', '--pdbqt_wildcard', help='Vina PDBQT poses filename.')
@click.option('-g', '--mol2_wildcard', help='Gold MOL2 results filenames.')
@click.option('-o', '--output', 'out_png', default='roc.png', help='Output PNG filename.')
def main(pdbqt_wildcard, mol2_wildcard,  active_regex, out_png):
    """Plots a Receiver-Operating-Characteristic curve for AutoDock Vina and
    GOLD virtual screening campaings.
    
    If the ACTIVE_REGEX match a Vina filename or Gold molecule name, it will be
    considered active.
    """
    try:
        active_regex = re.compile(active_regex)
    except:
        click.echo(click.style('Failed to compile ACTIVE_REGEX.', fg='red'))
        sys.exit(1)
    engines = []
    if pdbqt_wildcard:
        filenames = glob(pdbqt_wildcard)
        scores, actives = analyze_vina(filenames, active_regex)
        engines.append(('Vina', scores, actives))
    if mol2_wildcard:
        filenames = glob(mol2_wildcard)
        scores, actives = analyze_gold(filenames, active_regex)
        engines.append(('Gold', scores, actives))
    
    for (title, scores, actives) in engines:
        if len(scores) <= 1:
            click.echo(click.style(f"The wildcard must yield more than 1 molecule.", fg="red"))
            sys.exit(2)
        fpr, tpr, _ = roc_curve(actives, scores)
        auc = round(roc_auc_score(actives, scores), 2)

        plt.plot(fpr, tpr, label=f'{title}, auc={auc}')
        plt.legend(loc=4)
    plt.savefig(out_png)


if __name__ == '__main__':
    main()
