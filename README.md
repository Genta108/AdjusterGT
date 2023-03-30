# AdjusterGT

## parameter
### Switch of complement types
- JUMP_SW: 0or1
- STAY_SW: 0or1
- PARTS_SW: 0or1
- ORIGIN_SW: 0or1

### Result display
- DISPLAY: 0or1

### Variables for complement
- PARTS(8): Number of parts
- COLUMN:
- FRAME_RATE:
- JUMP_FIX_MAX:
- PINCH_FIX_MAX:
- PARTS_FIX_MAX:
- ROLL_FIX_MAX:
- FRAMEMAX:
- JUPMDET_WD:
- JUMPMOD_WD:
- STAYMOD_WD:
- STAYMOD_WD_MIN:
- PARTSMOD_WD:
- LIKELI_UP: 
- LIKELI_DOWN:
- PRECISION:
- VCLIM:
- DRLIM: 

### command example
- ./adjuster -group hippo -term 1 -f "filename" -pid A4 -cycle 1000 -sp 0
