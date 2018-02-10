CC=     g++
OPT=    -O3

all:
	$(CC) $(OPT) bgl_to_ped.cpp -o bgl_to_ped
	$(CC) $(OPT) ped_to_bgl.cpp -o ped_to_bgl
	$(CC) $(OPT) search.cpp -o search
	$(CC) $(OPT) ped_to_hom.cpp -o ped_to_hom
	$(CC) $(OPT) hom_to_ped.cpp -o hom_to_ped
	mv bgl_to_ped ped_to_bgl ped_to_hom hom_to_ped search ../../bin
clean:
	rm -f bgl_to_ped ped_to_bgl ped_to_hom hom_to_ped search
