SRC_MISC	=	cam_matrix.c	\
			cam_vector.c	\
			cam_list.c	\
			cam_3d_points_loaders.c	\
			cam_write_points_to_pgm.c

SRC_VIEWER	=	$(SRC_MISC)	\
			cam_3d_viewer.c

SRC_SFM		=	$(SRC_MISC)	\
			cam_project_2d_to_3d.c	\
			cam_project_3d_to_2d.c

OBJ_MISC	=	$(SRC_MISC:.c=.o)
OBJ_VIEWER	=	$(SRC_VIEWER:.c=.o)
OBJ_SFM		=	$(SRC_SFM:.c=.o)

NAME_VIEWER	=	viewer
NAME_SFM	=	sfm

CFLAGS		=	-W -Wall -ansi -pedantic -g3

CC		=	gcc

LIBS_VIEWER	=	-lGLU -lglut -lGL
LIBS_SFM	=	-lm

all: viewer sfm

viewer: $(OBJ_VIEWER)
	$(CC) $(OBJ_VIEWER) -o $(NAME_VIEWER) $(LIBS_VIEWER)

sfm: $(OBJ_SFM)
	$(CC) $(OBJ_SFM) -o $(NAME_SFM) $(LIBS_SFM)

clean:
	rm -f $(OBJ_VIEWER)
	rm -f $(OBJ_SFM)

fclean:	clean
	rm -f *~
	rm -f \#*\#
	rm -f $(NAME_VIEWER)
	rm -f $(NAME_SFM)

re : clean all