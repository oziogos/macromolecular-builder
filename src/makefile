

#################### Program ###################################


EXEC=builder3        

CFLAGS = -Wall -O3

SRC_DIR = .

IFLAGS = -I/$(SRC_DIR)

OBJS +=  \
	$(SRC_DIR)/builder2_euler.o \
	$(SRC_DIR)/builder2_initialize.o \
	$(SRC_DIR)/builder3.o \
	$(SRC_DIR)/builder_core_funct.o \
	$(SRC_DIR)/builder_loop_preamble.o \
	$(SRC_DIR)/contiguous_atom.o \
	$(SRC_DIR)/contiguous_bond.o \
	$(SRC_DIR)/find_unique.o \
	$(SRC_DIR)/find_unique_d.o \
	$(SRC_DIR)/grow_ch2.o \
	$(SRC_DIR)/grow_ch3.o \
	$(SRC_DIR)/grow_oh.o \
	$(SRC_DIR)/growth_reg_init.o \
	$(SRC_DIR)/lb.o \
	$(SRC_DIR)/permute.o \
	$(SRC_DIR)/resize_master_registries.o \
	$(SRC_DIR)/retrieve_mol_type_data.o \
	$(SRC_DIR)/rotate.o \
	$(SRC_DIR)/segments.o \
	$(SRC_DIR)/SYBYLtoUFF.o \
	$(SRC_DIR)/tokenize.o \
	$(SRC_DIR)/topo.o \
	$(SRC_DIR)/UFF_topo.o \
	$(SRC_DIR)/UFF_write_lammps.o \
	$(SRC_DIR)/write_xyz_preview.o \
	$(SRC_DIR)/topo_molecule_alloc_populate.o \
	$(SRC_DIR)/topo_general_backup.o \
	$(SRC_DIR)/alter_molecular_topo.o \
	$(SRC_DIR)/compute_local.o \
	$(SRC_DIR)/alter_general_topo.o \
	$(SRC_DIR)/tmap_topo_boundaries_init.o \
	$(SRC_DIR)/update_master.o \
	$(SRC_DIR)/topo_diagnostic.o \
	$(SRC_DIR)/remove_deleted.o \
	$(SRC_DIR)/grow_nh2.o \
	$(SRC_DIR)/initialize_tacticity.o \
	$(SRC_DIR)/apply_backbone_dihedral.o \
	$(SRC_DIR)/force_tacticity.o \
	$(SRC_DIR)/free_final.o \
	$(SRC_DIR)/free_locals.o \
	$(SRC_DIR)/free_molecules_loop.o \
	$(SRC_DIR)/minimize.o \
	$(SRC_DIR)/write_mol2_file.o \
	$(SRC_DIR)/calculate_tacticity.o \
	$(SRC_DIR)/update_neighbors_registry.o \
	$(SRC_DIR)/post_exec.o \
	$(SRC_DIR)/topo_neighbors.o \
	$(SRC_DIR)/linked_list_cell_initialize.o \
	$(SRC_DIR)/lj_pe_brute_force.o \
	$(SRC_DIR)/calc_head_list.o \
	$(SRC_DIR)/lj_pe_linked_list_cell.o \
	$(SRC_DIR)/LJ_params.o \
	$(SRC_DIR)/calc_wrapped_coords_n_arrays.o \
	$(SRC_DIR)/CBMCG.o \
	$(SRC_DIR)/pe_calc_module_v1.o \
	$(SRC_DIR)/calculate_dihedrals.o \
	$(SRC_DIR)/SYBYLtoDRE.o \
	$(SRC_DIR)/DRE_topo.o \
	$(SRC_DIR)/DRE_write_lammps.o \
   
STD_LIBS = -lm 
#################### Target #####################################

all: $(EXEC)

$(EXEC): $(OBJS)
	gcc $(CFLAGS) $(IFLAGS) $(OBJS) -o $@ $(STD_LIBS)	

clean: 
	rm -f $(OBJS) $(EXEC)

