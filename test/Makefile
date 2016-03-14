NAME	= viewer

CC	= g++

RM	= rm -f

SRCDIR	= src/

INCDIR	= inc/

CHLIBINC	= ../ConvexHull/inc/
QGDBSCANLIBINC =  ../inc/

OBJDIR	= obj/

CHLIBDIR 	= ../ConvexHull/build/
QGDBSCANLIBDIR = ../build/

CHLIBNAME = ConvexHull
QGDBSCANLIBNAME = QGDBSCAN

SRCS	= main.cpp \
	  Render.cpp
			
OBJS	= $(addprefix $(OBJDIR), $(SRCS:.cpp=.o))

INCCOMP	= -I $(INCDIR) -I $(CHLIBINC) -I $(QGDBSCANLIBINC)

CXXFLAGS	= -O2 -march=native -std=c++14

LIB	= -L$(CHLIBDIR) -L$(QGDBSCANLIBDIR) -l$(CHLIBNAME) -l$(QGDBSCANLIBNAME)

LDFLAGS	= -lsfml-graphics -lsfml-window -lsfml-system

dummy	:= 	$(shell test -d $(OBJDIR) || mkdir $(OBJDIR))

$(OBJDIR)%.o:	$(SRCDIR)%.cpp
	$(CC) $(CXXFLAGS) $(INCCOMP) -o $@ -c $<

$(NAME): $(OBJS)
	$(CC) -o $(NAME) $(OBJS) $(LDFLAGS) $(LIB)

all: $(NAME)

clean:
	$(RM) $(OBJS)

fclean:	clean
	$(RM) $(NAME)

re:	fclean all

.PHONY: all clean fclean re
