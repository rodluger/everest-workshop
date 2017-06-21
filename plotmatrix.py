import matplotlib.pyplot as pl
import numpy as np

def PlotMatrix(model, w, X):
    fig = pl.figure(figsize = (12,8))
    fig.subplots_adjust(wspace = 1)
    axm = pl.subplot2grid((1, 5), (0, 0), colspan = 1, rowspan = 1)
    axX = pl.subplot2grid((1, 5), (0, 1), colspan = 3, rowspan = 1)
    axw = pl.subplot2grid((1, 5), (0, 4), colspan = 1, rowspan = 1)
    axm.plot(model, np.arange(len(model)), color = 'r')
    axm.axis('off')
    axX.imshow(X, aspect = 'auto');
    axX.axis('off');
    axw.imshow(w, aspect = 0.2);
    axw.axis('off');
    fig.text(0.225,0.5,'=', fontsize = 50)
    fig.text(0.76,0.525,'.', fontsize = 100)