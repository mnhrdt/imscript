
struct point {
	int a, b;
} points[] = {
	{1,2},
	{2,3},
	{3,4}
};

int f()
{
	struct inner_point {
		int a, b;
		//int *s;
		struct stuff { int x, y, z; } *t;
	} the_inner_point = {
		1, 2,
		.t = (struct stuff []){
			{1,2,3},
			{1,2,3},
			{1,2,3}
			   }
	};

	//int *s = (int []){10, 20, 30, 40};
	//struct point *z = (struct point []){{1,2},{2,3},{4,5}};

	return f(),f();
}

int main(){return f();}
