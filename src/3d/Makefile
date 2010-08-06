SRC_MISC	=	cam_matrix.c		\
			cam_vector.c		\
			cam_list.c		\
			cam_3d_points_loaders.c	\
			cam_projection_matrix.c \
			cam_points_to_pgm.c	\
			cam_pgm_to_points.c	\
			cam_pgm_to_matrix.c

SRC_VIEWER	=	$(SRC_MISC)	\
			cam_3d_viewer.c

SRC_SFM		=	$(SRC_MISC)		\
			cam_epipoles.c		\
			cam_p_from_f.c		\
			cam_project_2d_to_3d.c	\
			cam_project_3d_to_2d.c	\
			sfm.c

SRC_NON_REGRESSION	=	$(SRC_MISC)		\
				cam_epipoles.c		\
				cam_p_from_f.c		\
				cam_project_2d_to_3d.c	\
				cam_project_3d_to_2d.c	\
				non_regression.c

SRC_EVALUATE_TRACKING	=	$(SRC_MISC)		\
				cam_evaluate_tracking.c

SRC_IMAGES_TRANSFORMATION =	$(SRC_MISC)				\
				cam_generate_images_transformations.c


OBJ_MISC			=	$(SRC_MISC:.c=.o)
OBJ_VIEWER			=	$(SRC_VIEWER:.c=.o)
OBJ_SFM				=	$(SRC_SFM:.c=.o)
OBJ_NON_REGRESSION		=	$(SRC_NON_REGRESSION:.c=.o)
OBJ_EVALUATE_TRACKING		=	$(SRC_EVALUATE_TRACKING:.c=.o)
OBJ_IMAGES_TRANSFORMATION	=	$(SRC_IMAGES_TRANSFORMATION:.c=.o)

NAME_VIEWER			=	viewer
NAME_SFM			=	sfm
NAME_NON_REGRESSION		=	non_regression
NAME_EVALUATE_TRACKING		=	evaluate_tracking
NAME_IMAGES_TRANSFORMATION	=	images_transformation

CFLAGS		=	-W -Wall -pedantic -ansi -Werror -g3

CC		=	gcc

LIBS_VIEWER			=	-lGLU -lglut -lGL
LIBS_SFM			=	-lm
LIBS_NON_REGRESSION		=	-lm
LIBS_EVALUATE_TRACKING		=	-lm
LIBS_IMAGES_TRANSFORMATION	=	-lm

all: viewer sfm non_regression evaluate_tracking images_transformation

viewer: $(OBJ_VIEWER)
	$(CC) $(OBJ_VIEWER) -o $(NAME_VIEWER) $(LIBS_VIEWER)

sfm: $(OBJ_SFM)
	$(CC) $(OBJ_SFM) -o $(NAME_SFM) $(LIBS_SFM)

non_regression:	$(OBJ_NON_REGRESSION)
	$(CC) $(OBJ_NON_REGRESSION) -o $(NAME_NON_REGRESSION) $(LIBS_NON_REGRESSION)

evaluate_tracking : $(OBJ_EVALUATE_TRACKING)
	$(CC) $(OBJ_EVALUATE_TRACKING) -o $(NAME_EVALUATE_TRACKING) $(LIBS_EVALUATE_TRACKING)

images_transformation : $(OBJ_IMAGES_TRANSFORMATION)
	$(CC) $(OBJ_IMAGES_TRANSFORMATION) -o $(NAME_IMAGES_TRANSFORMATION) $(LIBS_IMAGES_TRANSFORMATION)

clean:
	rm -f *~
	rm -f \#*\#
	rm -f $(OBJ_VIEWER)
	rm -f $(OBJ_SFM)
	rm -f $(OBJ_NON_REGRESSION)
	rm -f $(OBJ_EVALUATE_TRACKING)
	rm -f $(OBJ_IMAGES_TRANSFORMATION)

fclean:	clean
	rm -f $(NAME_VIEWER)
	rm -f $(NAME_SFM)
	rm -f $(NAME_EVALUATE_TRACKING)
	rm -f $(NAME_IMAGES_TRANSFORMATION)

re : clean all