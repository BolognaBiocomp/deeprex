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
