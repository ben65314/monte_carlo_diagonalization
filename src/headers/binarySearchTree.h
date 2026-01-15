#include "utilities.h"
#ifndef __binarySearchTree_h__
#define __binarySearchTree_h__
struct Node {
	sType key;
	Node* left;
	Node* right;
	Node() {
		left = right = NULL;
	}
	Node(sType item) {
		key = item;
		left = right = NULL;
	}

};

Node* insert(Node* root, Node* key, bool* alreadyThere);
const Node* search(const Node* root, uint64_t key);
void inorder(const Node* root);
void treeToVec(const Node* root, std::vector<uint64_t>* vec);
uint64_t totalNodes(Node* root);
Node* delNode(Node* root, uint64_t x);
int height(Node *root);

void balanceBST(std::vector<Node>* arr);
void printBT(const Node* node);

#endif
