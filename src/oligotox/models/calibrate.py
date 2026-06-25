from sklearn.calibration import CalibratedClassifierCV

def build_calibrator(estimator, method, cv):
    return CalibratedClassifierCV(estimator=estimator, method=method, cv=cv)
