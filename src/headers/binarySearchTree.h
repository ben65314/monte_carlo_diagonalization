#include "utilities.h"
#ifndef __binarySearchTree_h__
#define __binarySearchTree_h__
struct Node {
	uint64_t key;
	Node* left;
	Node* right;
	Node() {
		left = right = NULL;
	}
	Node(uint64_t item) {
		key = item;
		left = right = NULL;
	}

};

Node* insert(Node* root, uint64_t key, bool* alreadyThere);
Node* insert(Node* root, Node* key, bool* alreadyThere);
const Node* search(const Node* root, uint64_t key);
void inorder(const Node* root);
void treeToVec(const Node* root, std::vector<uint64_t>* vec);
uint64_t totalNodes(Node* root);
Node* delNode(Node* root, uint64_t x);
#endif
