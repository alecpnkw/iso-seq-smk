#! /bin/bash

bsub -P acc_rosenb16a -q interactive -n 2 -R "rusage[mem=4000] span[hosts=1]" -W 12:00 -Is /bin/bash
