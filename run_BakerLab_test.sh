#! /bin/bash

nice -n 10 python3 training_HotProtein.py \
	--path_for_outputs "/stor/work/Ellington/ProteinMPNN/HotProtein/Models/$class/model" \
	--path_for_training_data "/stor/work/Ellington/ProteinMPNN/HotProtein/Models/training_data"\
	--previous_checkpoint "/stor/work/Ellington/ProteinMPNN/HotProtein/ProteinMPNN_Data/ProteinMPNN_retrain/model_weights/epoch_last.pt"\
	--num_epochs 200\
