#include "binarySearchTree.h"

Node* insert(Node* node, uint64_t key, bool* alreadyThere) {
	// If the tree is empty, return a new node
	//std::cout<<"INSERTING:\tNODE:"<<node->key<<"\tKEY:"<<key->key<<std::endl;
    if (node == NULL){
		//std::cout<<"NULL, creating new node:"<<key<<std::endl;
        return new Node(key);    
	}
    // If the key is already present in the tree,
    // return the node
    if (node->key == key){ 
		//std::cout<<"ALREADY EXISTS"<<std::endl;
		*alreadyThere = true;
        return node;
	}
    // Otherwise, recur down the tree/ If the key
    // to be inserted is greater than the node's key,
    // insert it in the right subtree
    if (node->key < key){ 
		//std::cout<<"GOING RIGHT"<<std::endl;
        node->right = insert(node->right, key,alreadyThere);
	}
    // If the key to be inserted is smaller than 
    // the node's key,insert it in the left subtree
    else { 
		//std::cout<<"GOING LEFT"<<std::endl;
		//std::cout<<"LEFT NODE:"<<node->left<<std::endl;
		//std::cout<<"NODE:"<<node->key<<"\tKEY:"<<key->key<<std::endl;
        node->left = insert(node->left, key,alreadyThere);
	}
    // Return the (unchanged) node pointer
    return node;
}

Node* insert(Node* node, Node* key, bool* alreadyThere) {
	// If the tree is empty, return a new node
	//std::cout<<"INSERTING:\tNODE:"<<node->key<<"\tKEY:"<<key->key<<std::endl;
    if (node == NULL){
		//std::cout<<"NULL, creating new node:"<<key<<std::endl;
        return key;    
	}
    // If the key is already present in the tree,
    // return the node
    if (node->key == key->key){ 
		//std::cout<<"ALREADY EXISTS"<<std::endl;
		*alreadyThere = true;
        return node;
	}
    // Otherwise, recur down the tree/ If the key
    // to be inserted is greater than the node's key,
    // insert it in the right subtree
    if (node->key < key->key){ 
		//std::cout<<"GOING RIGHT"<<std::endl;
        node->right = insert(node->right, key,alreadyThere);
	}
    // If the key to be inserted is smaller than 
    // the node's key,insert it in the left subtree
    else { 
		//std::cout<<"GOING LEFT"<<std::endl;
		//std::cout<<"LEFT NODE:"<<node->left<<std::endl;
		//std::cout<<"NODE:"<<node->key<<"\tKEY:"<<key->key<<std::endl;
        node->left = insert(node->left, key,alreadyThere);
	}
    // Return the (unchanged) node pointer
    return node;
}

// function to search a key in a BST
const Node* search(const Node* root, uint64_t key) {
  
    // Base Cases: root is null or key 
    // is present at root
    if (root == NULL || root->key == key)
        return root;

    // Key is greater than root's key
    if (root->key < key)
        return search(root->right, key);

    // Key is smaller than root's key
    return search(root->left, key);
}

// A utility function to do inorder tree traversal
void inorder(const Node* root) {
    if (root != NULL) {
        inorder(root->left);
		std:: cout << root->key << " ";
        inorder(root->right);
    }
}


void treeToVec(const Node* root, std::vector<uint64_t>* vec) {
	if (root != NULL) {
		treeToVec(root->left,vec);
		vec->push_back(root->key);
		treeToVec(root->right,vec);
	}
}

// Function to get the count of nodes
// in complete binary tree
uint64_t totalNodes(Node* root)
{
    if (root == NULL)
        return 0;
 
    int l = totalNodes(root->left);
    int r = totalNodes(root->right);
 
    return 1 + l + r;
}


// Note that it is not a generic inorder
// successor function. It mainly works
// when right child is not empty which is 
// the case wwe need in BST delete
Node* getSuccessor(Node* curr)
{
    curr = curr->right;
    while (curr != NULL && curr->left != NULL)
        curr = curr->left;
    return curr;
}

// This function deletes a given key x from
// the give BST and returns modified root of
// the BST (if it is modified)
Node* delNode(Node* root, uint64_t x)
{

    // Base case
    if (root == NULL)
        return root;

    // If key to be searched is in a subtree
    if (root->key > x)
        root->left = delNode(root->left, x);
    else if (root->key < x)
        root->right = delNode(root->right, x);

    // If root matches with the given key
    else {

        // Cases when root has 0 children
        // or only right child
        if (root->left == NULL) {
            Node* temp = root->right;
            delete root;
            return temp;
        }

        // When root has only left child
        if (root->right == NULL) {
            Node* temp = root->left;
            delete root;
            return temp;
        }

        // When both children are present
        Node* succ = getSuccessor(root);
        root->key = succ->key;
        root->right = delNode(root->right, succ->key);
    }
    return root;
}


