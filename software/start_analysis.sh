#!/bin/sh
cd /home/ubuntu/clinical_applications
/home/ubuntu/.local/bin/jupyter notebook --port 80 --ip 0.0.0.0 --config /home/ubuntu/clinical_applications/software/.jupyter/jupyter_notebook_config.py
