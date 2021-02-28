import numpy as np
from genebench.datatypes import GeneMethodResult
from genebench.datatypes import GeneDiffInput
from genebench.diffmethods.base.diffmethod import DiffMethod
from genebench.utils import Utils


class ChDirPy(DiffMethod):
    def setup(self, config):
        self.logger = Utils.get_logger("ChDirPy")

    def post_proces_results(self, genes, values):
        values = np.absolute(values)
        max_score = np.max(values)
        if max_score != 0.0:
            values = values/max_score
        return self.sort(genes, values)

    def run(self, input: GeneDiffInput) -> GeneMethodResult:
        A = np.array(input.control)
        B = np.array(input.perturbed)
        genes = np.array(input.genes)
        A, B, genes = self._throw_away_rows_without_variance(A, B, genes)
        genes, scores = self._chdir(A, B, genes)
        return GeneMethodResult.from_separate_lists(genes, scores)

    def _chdir(self, A, B, genes, r=1):
        X = np.concatenate((A, B), axis=1).T
        (rowCount, colCount) = np.shape(X)
        if 20 > rowCount-1:
            maxComponentsNum = rowCount - 1
        else:
            maxComponentsNum = 20
        scores, loadings, explained_var = self.nipals(X,
                                                      maxComponentsNum,
                                                      1e5,
                                                      1e-4)
        scores = scores.T
        loadings = loadings.T

        captured_variance = 0
        for i in range(len(explained_var)):
            captured_variance += explained_var[i]
            if captured_variance > 0.95:
                break
        scores = scores[0: i+1]
        loadings = loadings[0: i+1]
        scores = scores.T
        loadings = loadings.T
        meanvec = np.mean(B, axis=1) - np.mean(A, axis=1)
        Dd = np.dot(scores.T, scores)/rowCount
        Dd = np.diag(np.diag(Dd))
        sigma = np.mean(np.diag(Dd))
        shrunkMats = np.linalg.inv(np.dot(r, Dd) + sigma*(1-r)*np.eye(np.shape(Dd)[0]))
        b = np.dot(loadings, np.dot(shrunkMats, np.dot(loadings.T, meanvec)))
        b /= np.linalg.norm(b)
        return (genes, b.tolist())

    def nipals(self, X, a, it=100, tol=1e-4):
        X = np.array(X)
        (obsCount, varCount) = np.shape(X)
        Xh = X - np.tile(np.mean(X, axis=0), (obsCount, 1))
        T = np.zeros((obsCount, a))
        P = np.zeros((varCount, a))
        pcvar = np.zeros(varCount)
        varTotal = np.sum(np.var(Xh, axis=0))
        currVar = varTotal
        nr = 0

        for h in range(a):
            th = np.reshape(Xh[:, 0], (obsCount, -1))
            ende = False

            while not ende:
                nr = nr + 1
                ph = np.dot(Xh.T, th)/np.dot(th.T, th)
                ph = ph/np.sqrt(np.dot(ph.T, ph))
                thnew = np.dot(Xh, ph)/np.dot(ph.T, ph)
                prec = np.dot((thnew-th).T, (thnew-th))
                th = thnew

                if prec <= np.power(tol, 2):
                    ende = True
                if it <= nr:
                    ende = True
                    self.logger.warning('Iteration stops without convergence')

            Xh = Xh - np.dot(th, ph.T)
            T[:, h] = th[:, 0]
            P[:, h] = ph[:, 0]
            oldVar = currVar
            currVar = np.sum(np.var(Xh, axis=0))
            pcvar[h] = (oldVar - currVar)/varTotal
            nr = 0

        return T, P, pcvar

    def _throw_away_rows_without_variance(self, A, B, genes):
        """Discards rows without variance.
        """
        X = np.concatenate((A, B), axis=1)
        keep_rows = np.std(X, axis=1).nonzero()[0]
        return A[keep_rows, ], B[keep_rows, ], genes[keep_rows]
