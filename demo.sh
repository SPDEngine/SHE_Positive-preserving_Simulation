#!/usr/bin/env bash
# #!/bin/bash

she-sim --config ./she-sim-config.yaml
ffmpeg -i Demo_2-1_2d.mp4 Demo_2-1_2d.gif
ffmpeg -i Demo_2-2_2d.mp4 Demo_2-2_2d.gif
