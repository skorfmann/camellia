SRC_MISC	=	cam_matrix.c	\
			cam_list.c

SRC_VIEWER	=	$(SRC_MISC)	\
			cam_3d_viewer.c

SRC_SFT		=	$(SRC_MISC)	\
			cam_project_3d_to_2d.c

OBJ_MISC	=	$(SRC_MISC:.c=.o)
OBJ_VIEWER	=	$(SRC_VIEWER:.c=.o)
OBJ_SFT		=	$(SRC_SFT:.c=.o)

NAME_VIEWER	=	viewer
NAME_SFT	=	sft

CFLAGS		=	-W -Wall -ansi -pedantic -g3

CC		=	gcc

LIBS_VIEWER	=	-lGLU -lglut -lGL

all: viewer sft

viewer: $(OBJ_VIEWER)
	$(CC) $(OBJ_VIEWER) -o $(NAME_VIEWER) $(LIBS_VIEWER)

sft: $(OBJ_SFT)
	$(CC) $(OBJ_SFT) -o $(NAME_SFT)

clean:
	rm -f $(OBJ_VIEWER)
	rm -f $(OBJ_SFT)

fclean:	clean
	rm -f *~
	rm -f \#*\#
	rm -f $(NAME_VIEWER)
	rm -f $(NAME_SFT)

re : clean all