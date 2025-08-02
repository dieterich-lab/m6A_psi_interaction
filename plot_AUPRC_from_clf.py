import os
from glob import glob
import pickle
import matplotlib as mpl
cm = 1/2.54  # centimeters in inches
gr = 1.618
dpi = 1200
mpl.rcParams['figure.dpi'] = dpi
mpl.rcParams['savefig.dpi'] = dpi
mpl.rcParams['font.size'] = 8
mpl.rcParams['legend.fontsize'] = 6
mpl.rcParams['xtick.labelsize'] = 8
mpl.rcParams['ytick.labelsize'] = 8
mpl.rcParams['xtick.major.size'] = 4
mpl.rcParams['ytick.major.size'] = 4
mpl.rcParams['lines.linewidth'] = 1
# mpl.rcParams['font.family'] = 'Arial'
FMT = 'svg'
fig_kwargs = dict(format=FMT, bbox_inches='tight', dpi=dpi, transparent=True)
import matplotlib.pyplot as plt


clf_dir = "/home/achan/prj/TRR319_RMaP_BaseCalling/Adrian/classifiers/m6A"
clf_paths = glob(os.path.join(clf_dir, "*.pkl"))

out_dir = "/home/achan/prj/TRR319_RMaP_BaseCalling/Adrian/AUC/m6A"
os.makedirs(out_dir, exist_ok=True)

for this_clf_path in clf_paths:
    this_motif = os.path.basename(this_clf_path).rstrip('.pkl')
    out_img_path = os.path.join(out_dir, f"{this_motif}.{FMT}")

    with open(this_clf_path, "rb") as handle_clf:
        clf = pickle.load(handle_clf)


    plt.figure(figsize=(5*cm, 5*cm))
    plt.plot(clf.recall, clf.precision, '-')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.xlim([0, 1.05])
    plt.ylim([0, 1.05])
    plt.title(f'{clf.motif}\n{len(clf.train_nts)} train NTs, {len(clf.test_nts)} test NTs\n AUC = {clf.auc:.2f}')
    plt.savefig(out_img_path, **fig_kwargs)
    plt.close('all')