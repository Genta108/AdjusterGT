# AdjusterGT

## parameter
### Switch of complement types
- JUMP_SW: 0or1
- STAY_SW: 0or1
- PARTS_SW: 0or1
- ORIGIN_SW: 0or1

### Result displaying
- DISPLAY: 0or1

### Variables for complement (default value)
- PARTS(8): the number of rat's parts
- COLUMN(25): columns of result matrix
- FRAME_RATE(15): video fps
- JUMP_FIX_MAX(10): times of modification for jumping plot in a cycle
- PINCH_FIX_MAX(100): times of pinching modification for long term jump in a cycle
- PARTS_FIX_MAX(1): times of modification for outlier parts in a cycle
- ROLL_FIX_MAX(10): times of rolling modification for long term jump in a cycle
- FRAMEMAX(486001): the number of video frames
- JUPMDET_WD(FRAME_RATE x 3): the number of records in the window for the detection of jumping plot
- JUMPMOD_WD(FRAME_RATE): the number of records in the window for the modification of jumping plot
- STAYMOD_WD(FRAME_RATE x 3): the number of records in the window for the modification of long term jump
- STAYMOD_WD_MIN(FRAME_RATE x 100): the minimal number of 
- PARTSMOD_WD(FRAME_RATE):
- LIKELI_UP(0.1): 
- LIKELI_DOWN(1):
- PRECISION(1):
- VCLIM(30):
- DRLIM(300): 

### command example
- ./adjuster -group hippo -term 1 -f "filename" -pid A4 -cycle 1000 -sp 0
