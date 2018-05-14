import numpy as np
from matplotlib.colors import Normalize
import matplotlib.cm as cm


class ColourScale(object):
    """ColourScale for plotting clones and mutations"""

    def __init__(self, css, by_clone_type=True, by_ns=True, all_clones_noisy=False, by_initial=False):
        self.by_clone_type = by_clone_type
        self.by_ns = by_ns
        self.all_noise = all_clones_noisy
        self.by_initial = by_initial
        self.css = css  # colourmap if only one for all clones, otherwise a dictionary

    def get_colour(self, rate, clone_type, ns, initial):
        if self.by_clone_type:
            if self.by_ns:
                if self.by_initial:
                    cs = self.css[(clone_type, ns, initial)]
                else:
                    cs = self.css[(clone_type, ns)]
            else:
                if self.by_initial:
                    cs = self.css[(clone_type, initial)]
                else:
                    cs = self.css[clone_type]
        elif self.by_ns:
            if self.by_initial:
                cs = self.css[(ns, initial)]
            else:
                cs = self.css[ns]
        elif self.by_initial:
            cs = self.css[initial]
        else:
            cs = self.css

        if (not ns and not initial) or self.all_noise:  # Add a bit of randomness to see synonymous mutations
            noise = np.random.normal(loc=np.random.choice([0.97, 1.03]), scale=0.01)
        else:
            noise = 0
        return cs(rate + noise)


COLOURSCALE_BAL = ColourScale(
    by_clone_type=False, by_ns=True, all_clones_noisy=True,
    css={True: cm.Greens,
         False: cm.ScalarMappable(norm=Normalize(vmin=0, vmax=2), cmap=cm.YlOrBr).to_rgba}
)

COLOURSCALE_MUT = ColourScale(
    by_clone_type=False, by_ns=True, all_clones_noisy=False,
    css={True: cm.Greens,
         False: cm.ScalarMappable(norm=Normalize(vmin=0, vmax=2), cmap=cm.YlOrBr).to_rgba}
)


def random_colour(rate):
    return cm.gnuplot(np.random.random())


def get_colourscale_with_random_mutation_colour(max_fitness):
    diff_ = max_fitness - 1
    cs = ColourScale(
        by_clone_type=True, by_ns=False, all_clones_noisy=False, by_initial=True,
        css={
            (1, True): cm.Greens,
            (1, False): cm.ScalarMappable(norm=Normalize(vmin=1 - 5*diff_, vmax=max_fitness), cmap=cm.Greens).to_rgba,
            (0, True): cm.ScalarMappable(norm=Normalize(vmin=0, vmax=3),
                                      cmap=cm.YlOrBr).to_rgba,
            (0, False): random_colour
            }
    )
    return cs