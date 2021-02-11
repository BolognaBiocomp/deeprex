import os
try:
    DEEPREX_ROOT = os.environ['DEEPREX_ROOT']
except:
    logging.error("DEEPREX_ROOT environment varible is not set")
    logging.error("Please, set and export DEEPREX_ROOT to point to deeprex root folder:")
    logging.error("$ export DEEPREX_ROOT=/path/to/deeprex")
    sys.exit(1)

HHBLITS_ITERATIONS = 2

DEEPREX_MODEL_FILE = os.path.join(DEEPREX_ROOT, "data", "model_final.h5")

KD = {'A': 1.8, 'C': 2.5, 'E': -3.5, 'D': -3.5,
      'G': -0.4, 'F': 2.8, 'I': 4.5, 'H': -3.2,
      'K': -3.9, 'M': 1.9, 'L': 3.8, 'N': -3.5,
      'Q': -3.5, 'P': -1.6, 'S': -0.8, 'R': -4.5,
      'T': -0.7, 'W': -0.9, 'V': 4.2, 'Y': -1.3, 'X': 0.0}
