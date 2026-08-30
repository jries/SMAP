# CLI

If installed as a package, the console script is `comet`. Otherwise:
```
python cli_interface.py   --input in.csv --output out.csv --format csv   --segmentation_mode 2 --segmentation_var 60   --initial_sigma_nm 100 --target_sigma_nm 1 --boxcar_width 1
```

**Important flags**
- `--segmentation_mode {0,1,2}`: windows / locs per window / frames per window
- `--segmentation_var`: depends on mode (e.g., frames per window for mode 2)
- `--initial_sigma_nm`, `--target_sigma_nm`
- `--boxcar_width`
- `--interpolation_method {cubic,catmull-rom}`

> Code default `initial_sigma_nm=100` differs from older README values.