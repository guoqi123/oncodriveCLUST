# Analytical p_value calculation

# Import modules
import numpy as np
from scipy.stats import gaussian_kde
from sklearn.neighbors import KernelDensity
from sklearn.model_selection import GridSearchCV


class AnalyticalPvalue:
    """Class to calculate the analytical pvalue of elements"""

    def __init__(self, jobs=1):
        """Initialize the AnalyticalPvalue class
        :param jobs: int, cores to use in the computation
        """
        self.jobs = jobs

    def _get_best_estimator(self, expected):
        """Calculate the best bandwidth estimator.
        :param expected: list of floats, list of expected values
        :return: float, best bandwidth estimator
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

    def _get_min_analytical_pvalue(self, gkde, up, down, calls=0, pvals=None):
        """Get the best pvalue with a resolution > 0.
        :param gkde: gaussian_kde object
        :param up: numeric, higher score
        :param down: numeric, lower score
        :param calls: int, number of times the function has been called
        :param pvals: list of pvalues
        :return: float, pvalue > 0
        """
        if pvals is None:
            pvals = []
        calls += 1
        mid = down + ((up - down) / 2)
        analytical_pvalue = gkde.integrate_box_1d(mid, float('inf'))
        if analytical_pvalue > 0:
            pvals.append(analytical_pvalue)
        if calls == 100:
            return pvals[-1]
        if analytical_pvalue > 0:
            return self._get_min_analytical_pvalue(gkde, up=up, down=mid, calls=calls, pvals=pvals)
        else:
            return self._get_min_analytical_pvalue(gkde, up=mid, down=down, calls=calls, pvals=pvals)

    def _calculate_analytical_pvalue(self, observed, expected):
        """Calculate the analytical pvalue of an element
        :return: float, analytical pvalue
        """
        gkde = gaussian_kde(expected)
        bandwidth = self._get_best_estimator(expected)
        gkde.set_bandwidth(bw_method=bandwidth)
        analytical_pvalue = gkde.integrate_box_1d(observed, float('inf'))
        if analytical_pvalue == 0:
            analytical_pvalue = self._get_min_analytical_pvalue(gkde, up=observed, down=0)
        return analytical_pvalue, bandwidth

    def calculate(self, observed, expected):
        """Calculate the analytical pvalue
        :param observed: float, observed value
        :param expected: list of floats, list of expected values
        :return: float, analytical pvalue
        """
        expected = np.array(expected)
        analytical_pvalue, bandwidth = self._calculate_analytical_pvalue(observed, expected)
        return analytical_pvalue, bandwidth