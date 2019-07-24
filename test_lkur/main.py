
from lightkurve import KeplerTargetPixelFile
tpf = KeplerTargetPixelFile.from_archive("Kepler-10", quarter=3, quality_bitmask='hardest')

