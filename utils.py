import numpy as np
from sklearn.metrics import mean_squared_error, r2_score
import seaborn as sns
import matplotlib.pyplot as plt



def start_plot(figsize=(10, 8), style = 'whitegrid', dpi=100):
    fig = plt.figure(figsize=figsize, dpi=dpi)
    gs = fig.add_gridspec(1,1)
    plt.tight_layout()
    with sns.axes_style(style):
        ax = fig.add_subplot(gs[0,0])
    return ax

def R2_plot(reconstruct_y_pred, reconstruct_y_test, error, accuracy, prop_name_T, ax=None):
    if ax is None:
        ax = start_plot(style='darkgrid', dpi=180)
    
    sns.set(rc={'font.family': 'Times New Roman'})
    ax.scatter(reconstruct_y_test.reshape(-1), reconstruct_y_pred.reshape(-1), color='darkorange', edgecolor='navy',
               label=r'$R^2$' + ": %.4f" % r2_score(y_true=reconstruct_y_test, y_pred=reconstruct_y_pred) + '\n' +
                     "MAPE" + ": %.4f" % error + '\n' +
                     "Accuracy" + ": %.4f" % accuracy)
    ymin = min(np.min(reconstruct_y_test), np.min(reconstruct_y_pred)) - 0.1
    ymax = max(np.max(reconstruct_y_test), np.max(reconstruct_y_pred)) + 0.1
    lim = [ymin, ymax]
    ax.set_xlim(lim)
    ax.set_ylim(lim)
    ax.plot(lim, lim, c='brown', ls='--', label=r'$y=\hat y, $' + 'identity')
    ax.legend(loc='best', frameon=True, shadow=True, fontsize=15)
    plt.xlabel('TRUE %s' % prop_name_T, fontsize=20, font="Times New Roman")
    plt.ylabel('PREDICTED %s' % prop_name_T, fontsize=20, font="Times New Roman")
    plt.xticks(fontsize=20, font="Times New Roman")
    plt.yticks(fontsize=20, font="Times New Roman")
    plt.title('%s Testing: %d' % (prop_name_T, len(reconstruct_y_test)), fontsize=20, font="Times New Roman")
    return ax


def melting_points(name):
    if name == "EC":
        return 36.4
    elif name == "PC":
        return -48.8