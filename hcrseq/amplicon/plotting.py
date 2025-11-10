import seaborn as sns

def temporal_plot(df,x="timepoint",
                  y="value",
                  col="pathway",
                  hue="cell_line",
                  sd='sd',
                  units = None,
                  estimator = 'mean',
                  draw_individual_errorbars=True,
                  draw_individual_points = False,
                  col_wrap=3,
                  height=2.5,
                  **kwargs):

    if (units is not None):
        estimator = None

    g= sns.relplot(x=x,y=y,
                   col=col,
                   hue=hue,
                   kind ="line",
                   units = units,
                   estimator = estimator,
                   data=df,
                   col_wrap=col_wrap,height=height,
                   facet_kws={'sharey': False},**kwargs)

    if draw_individual_points:
        for ax, (_, data) in zip(g.axes.flat, g.facet_data()):
            if not data.empty:
                sns.scatterplot(
                    x=x,
                    y=y,
                    hue=hue,
                    data=data,
                    ax=ax,  # Plot on the current subplot
                    #dodge=True,  # Separate points by hue, like the lines
                    alpha=0.6,  # Make points slightly transparent
                    size=4,  # Adjust point size
                    legend=False  # Avoid adding a duplicate legend
                )

    # Add error bars
    if (units is not None) and draw_individual_errorbars:
        for ax,(ind,data) in zip(g.axes,g.facet_data()):
            for ind2,datag in data.groupby(hue):
                ax.errorbar(x=datag[x],y =datag[y],yerr=datag[sd],fmt='none',capsize=2,color='grey')

    for pathway in df['pathway'].unique():
        if (pathway not in ['HR','NHEJ','NER']) and (not pathway.startswith('MMEJ')):
            g.axes_dict[pathway].set_ylim((0,1))

    # Remove pathway= from labels
    g.set_titles(col_template="{col_name}", row_template="{row_name}")

    for i, ax in enumerate(g.axes.flat):
        # Only show y-label on the first column of the grid
        if i % col_wrap != 0:
            ax.set_ylabel('')
        else:
            ax.set_ylabel('Repair Capacity')

    return(g)