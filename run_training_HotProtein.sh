#! /bin/bash

s2c2=(class_0_1_2 class_3_4)
s2c5=(class_0 class_1 class_2 class_3 class_4)

for class in "${s2c2[@]}"; do
	cp /stor/work/Ellington/ProteinMPNN/HotProtein/Models/$class/*txt /stor/work/Ellington/ProteinMPNN/HotProtein/Models/training_data/
	cp /stor/work/Ellington/ProteinMPNN/HotProtein/Models/$class/list.csv /stor/work/Ellington/ProteinMPNN/HotProtein/Models/training_data/
	
	nice -n 10 python3 training_HotProtein.py \
		--path_for_outputs "/stor/work/Ellington/ProteinMPNN/HotProtein/Models/$class/model" \
		--path_for_training_data "/stor/work/Ellington/ProteinMPNN/HotProtein/Models/training_data"\
		--previous_checkpoint "/stor/work/Ellington/ProteinMPNN/HotProtein/ProteinMPNN_Data/ProteinMPNN_retrain/model_weights/epoch_last.pt"\
		--num_epochs 200\
		--save_model_every_n_epochs 5
done

