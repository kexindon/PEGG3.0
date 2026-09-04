# create the environment and install PEGG
conda create -n pegg_env python=3.9 -y && conda activate pegg_env && pip install pegg
# or update PEGG
# conda activate pegg_env && pip install --upgrade --upgrade-strategy only-if-needed pegg

# to use web tool to visualize pegRNA design
conda activate pegg_env
cd ~/Documents/GitHub/prime-editing-visualization

pip install -r requirements.txt # first time use only
python run.py

# open http://127.0.0.1:5050 in your browser