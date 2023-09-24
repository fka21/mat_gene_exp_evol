#!/bin/bash

orthofinder -t 15 -a 15 -M msa -S diamond_ultra_sens -f /home/s9/tsai/work/Ferenc/Ortho/primary_transcripts/
orthofinder -t 15 -a 15 -M msa -S diamond_ultra_sens -s /home/s9/tsai/work/Ferenc/Ortho/rotl_species.nwk -f /home/s9/tsai/work/Ferenc/Ortho/primary_transcripts/

