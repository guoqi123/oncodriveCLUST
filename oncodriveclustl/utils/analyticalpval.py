"""
Contains class AnalyticalPvalue for analytical p_value calculation
"""

# Import modules
import numpy as np
from scipy.stats import gaussian_kde
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV

import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)


class AnalyticalPvalue:
    """Class to calculate the analytical pvalue of elements"""

    def __init__(self, jobs=1):
        """Initialize the AnalyticalPvalue class

        Args:
            jobs (int): cores to use in the computation

        Returns:
            None
        """

        self.jobs = jobs
        self.gkde = None
        self.bandwidth = None

    def _get_best_estimator(self, expected):
        """Calculate the best bandwidth estimator.

        Args:
            expected: list of floats, list of expected values

        Returns:
            float, best bandwidth estimator
        """

        # parameters of the fitting
        params = {
            'bandwidth': np.logspace(-1, 0, 10),
            'kernel': ['gaussian']
        }
        grid = GridSearchCV(KernelDensity(), params, n_jobs=self.jobs)
        data_newaxis = expected[:, np.newaxis]
        grid.fit(data_newaxis)
        return grid.best_estimator_.bandwidth

    def _get_min_analytical_pvalue(self, up, down, calls=0, pvals=None):
        """Get the best pvalue with a resolution > 0.

        Args:
            up (numeric): higher score
            down (numeric): lower score
            calls (int): number of times the function has been called
            pvals (list): list of pvalues

        Returns:
            float, pvalue > 0
        """

        if pvals is None:
            pvals = []
        calls += 1
        mid = down + ((up - down) / 2)
        analytical_pvalue = self.gkde.integrate_box_1d(mid, float('inf'))
        if analytical_pvalue > 0:
            pvals.append(analytical_pvalue)
        if calls == 100:
            return pvals[-1]
        if analytical_pvalue > 0:
            return self._get_min_analytical_pvalue(up=up, down=mid, calls=calls, pvals=pvals)
        else:
            return self._get_min_analytical_pvalue(up=mid, down=down, calls=calls, pvals=pvals)

    def calculate_bandwidth(self, expected, pseudocount=1):
        """Calculate the best estimator
        Args:
            expected (iterable): array of simulated mutations
            pseudocount (int): value to add when expected are zeroes

        Returns:
            None

        """

        # Add a pseudocount to get rid of LinAlgError("singular matrix")
        if sum(expected) == 0:
            expected = np.array([0] * 999 + [pseudocount])
        else:
            expected = np.array(expected)
        try:
            self.gkde = gaussian_kde(expected)
        except np.linalg.linalg.LinAlgError as e:
            expected[-1] += 10**(-14)
            self.gkde = gaussian_kde(expected)
        self.bandwidth = self._get_best_estimator(expected)
        self.gkde.set_bandwidth(bw_method=self.bandwidth)

    def _calculate_analytical_pvalue(self, observed):
        """Calculate the analytical pvalue of an element

        Args:
           observed (float): observed value

        Returns:
            float, analytical pvalue

        """

        analytical_pvalue = self.gkde.integrate_box_1d(observed, float('inf'))
        if analytical_pvalue == 0:
            analytical_pvalue = self._get_min_analytical_pvalue(up=observed, down=0)
        return analytical_pvalue

    def calculate(self, observed):
        """Calculate the analytical pvalue

        Args:
           observed (float): observed value

        Returns:
            analytical_pvalue (float): analytical pvalue
        """

        analytical_pvalue = self._calculate_analytical_pvalue(observed)
        return analytical_pvalue
