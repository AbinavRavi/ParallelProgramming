#include "familytree.h"


void traverse_tree(tree *node)
{
	if(node != NULL)
	{
		node->IQ = compute_IQ(node->data);
		genius[node->id] = node->IQ;

		#pragma omp task
		{
			traverse_tree(node->right);
		}

		#pragma omp task
		{
			traverse_tree(node->left);
		}
	}
}

void traverse(tree *node, int numThreads)
{

	// put your code in here
	#pragma omp parallel num_threads(numThreads)
	{
		#pragma omp single
		{
		traverse_tree(node);
		}
	}
}
