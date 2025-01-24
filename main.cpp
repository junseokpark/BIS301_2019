#include <iostream>
#include <map>
#include <list>
#include <string>
#include <fstream>
#include <sstream>

#include <vector>
#include <math.h>

using namespace std;


/***********************************************************************************************
   Part 1. Structure for Decision Tree ( for pre-lab, you don't need to change this part )
***********************************************************************************************/
struct Data // to store microarray data sample and gene list
{
    list<unsigned int> sampleList; // sample list
    list<string> clsGeneList; // candidate gene list for classification
};

struct Class // to store classification information of samples
{
    string clsResult; // class of a node - normal or disease; ( inputted only when a node is leaf node, default is null;)
    string clsGene; // gene for best classification
    int cutoffValue; // cutoff value of the clsGene for classification
};

struct Node // to organize binary tree structure
{
    Data d;
    Class cls;
    Node* high; // child node (expression => cls_value);
    Node* low; // child node (expression < cls_value);
};

std::string trim(string& str, const string& whitespace = " \t")
{
    unsigned int strBegin = str.find_first_not_of(whitespace);
    if (strBegin == std::string::npos)
        return ""; // no content

    unsigned int strEnd = str.find_last_not_of(whitespace);
    unsigned int strRange = strEnd - strBegin + 1;

    return str.substr(strBegin, strRange);
};
// split string by delimeter into vector.
using string_vector = vector<string>;
string_vector split(string str, char delimiter) {
    vector<string> internal;
    stringstream ss(str); // Turn the string into a stream.
    string tok;

    while (getline(ss, tok, delimiter)) {
        //Remove whitespaces
        tok = trim(tok);
        internal.push_back(tok);
    }

    return internal;
};

double entropycal(int fOrM, int lAtT) {
    double FoRm = (double)fOrM; double LaTt = (double)lAtT;
    if (fOrM == 0 && lAtT == 0) return 0.00;
    else if (fOrM == 0) return 0.00;
    else if (lAtT == 0) return 0.00;
    else return -(FoRm / (FoRm + LaTt)) * log2(FoRm / (FoRm + LaTt)) -(LaTt / (FoRm + LaTt)) * log2(LaTt / (FoRm + LaTt));
};


/***********************************************************************************************
   Part 2. Declaration of class DecisionTree ( for pre-lab, you don't need to change this part )
***********************************************************************************************/

using string_list = list<string>;
class DecisionTree
{
private:
    map<unsigned int, bool> sampleCls; // Sample property: normal(TRUE) or cancer(FALSE)
    map<pair<unsigned int, string>, double> exp; // microarray expression data

public:
    Node dc; // will be used as a root node of a decision tree

    // You may complete following functions at implementation part (part 3);
    void readDataFromFile(string filename); // read & save microarray data into 'exp' variable of this class.
    string_list readDEGFromFile(string filename); // read & return differentially expressed gene list.
    double findCutoffValue(string gene, list<unsigned int> sampleList); // find a best cutoff value of classification for an input gene and input samples.
    double getInfoGain(string gene, double cutoffValue, list<unsigned int> sampleList); // calculate information gain for an input gene and sample list.
    void testReadFile(  );
    void testDEGList (list<string> deg);
    void testInfoGain(string gene, list<double> cutoffValueList, list<unsigned int> sampleList); // Checking information gain value for
    void testCutoff(string gene, list<unsigned int> sampleList); // Checking finding elements from lists
};
/* [HINT] for implementation of methods in class DecisionTree

 @ sampleCls variable
┌────────┬──────────┐
│Sample #│Sample cls│
├────────┼──────────┤
│        │          │
│  ...   │   ...    │
│        │          │
└────────┴──────────┘
 --> you can find sample class(normal/disease) for input, sample #, by calling pre-defined function find() of map container.


 @ exp variable
┌─────────────────────┬────────────────┐
│┌────────┬──────────┐│                │
││Sample #│Gene name ││Expression Value│
│└────────┴──────────┘│                │
├─────────────────────┼────────────────┤
│                     │                │
│          ...        │      ...       │
│                     │                │
└─────────────────────┴────────────────┘
--> you can find expression value for input gene of input sample # by find() function of map container with pair<> container.
*/


/***********************************************************************************************
   Part 3. Definition of methods(functions) in class `DecisionTree`
   ------ You should implement some of the methods(functions) in this part for Pre-lab.
***********************************************************************************************/
void DecisionTree::readDataFromFile(string filename){
    /* Goal of this method (function):
    // Reading expression data from file should be implemented.
    // 1. open the data file with `filename`
    // 2. read first line (sample classes)
    // 3. read following lines (expression values for a gene of each sample)
    // 4. insert data from file to proper lists
    //// The sample index and class of each sample should be saved at `sampleCls`.
    //// The expression values for genes of sample is saved at `exp`.
    // 5. return nothing (void)

    // TIP! : Each item of a line has been separated with `\t` (tab).
    // TIP! : To check your expression data, try calculate sum of expression value of each sample.
    */

    ifstream temp_file;
    string oneline; string onecell;
    temp_file.open( filename);
    if (!temp_file.is_open()) {
        cout << "Error : It can't be opend" << endl;
        exit(0);
    }
    getline(temp_file, oneline);
    istringstream ABC(oneline); getline(ABC, onecell, '\t');
    int cnt = 1;
    while (getline(ABC, onecell, '\t')) {
        if (onecell == "non-tumor") sampleCls.insert(pair<int, bool>(cnt, true));
        else sampleCls.insert(pair<int, bool>(cnt, false));
        cnt++;
    }
    while (getline(temp_file, oneline)) {
        string Name_p;
        istringstream ABC(oneline); getline(ABC, Name_p, '\t');
        int cnt = 1;
        while (getline(ABC, onecell, '\t')) {
            exp.insert(pair<pair<unsigned int, string>, double>(pair<unsigned int, string>(cnt, string(Name_p)), stod(onecell)));
            cnt++;
        }
    }


    testReadFile();

};
void DecisionTree::testReadFile(  ){
    int cnt;
    /* Test1
	// This is for testing your code.
	// Printing out index and class of it for each sample from `sampleCls` list */
    cnt = 0;
    map<unsigned int, bool>::iterator it = sampleCls.begin();
    cout << "TEST1: print list variable" << endl;
    for(it = sampleCls.begin(); it!=sampleCls.end(); it++){
        cout << "Sample " << (*it).first << ": " << ((*it).second ? "Normal" : "Disease")  << endl;
        if (++cnt >= 20){
            cout << "..." << endl;
            break;
        }
    }
    cout << endl;

    /* Test2
    // This is for testing your code.
    // Printing out expression value of cetain sample and gene	from exp list */
    cnt = 0;
    map<pair<unsigned int, string>, double>::iterator it2 = exp.begin();
    cout << "TEST2: print map variable" << endl;
    for(it2 = exp.begin(); it2!=exp.end(); it2++){
        cout << "(Sample " << (*it2).first.first << ", " << (*it2).first.second << "): " << (*it2).second  << endl;
        if (++cnt >= 20){
            cout << "..." << endl;
            break;
        }
    }
    cout << endl;
};

string_list DecisionTree::readDEGFromFile(string filename){
    /* Goal of this method (function):
    // Reading gene list from file should be implented and the list should be returned.
    // 1. open the data file with `filename`
    // 2. read lines (the name of DEGs)
    // 3. insert data from file to list of string. e.g. 'deg'
    // 4. return the list of name of DEG (string)

    // TIP! : Check the all of name of DEG is included in the name of genes of `exp` list
    */

    list<string> deg;

    ifstream temp_file;
    string oneline; string onecell;
    temp_file.open(filename);
    if (!temp_file.is_open()) {
        cout << "Error : It can't be opend" << endl;
        exit(0);
    }
    int cnt = 0;
    while (getline(temp_file, oneline)) {
        istringstream ABC(oneline); getline(ABC, onecell, '\t');
        if ( cnt++ >= 1)
            deg.push_back(string(trim(onecell) ));
    }

    return deg;
};
void DecisionTree::testDEGList (list<string> deg){
    int cnt;
    /* Test3
	// This is for testing your code.
	// Printing out all of the name of DEGs */
    cnt = 0;
    list<string>::iterator it;
    cout << "TEST3: print DEG list" << endl;
    for(it= deg.begin(); it!=deg.end(); it++){
        cout << *it << " | " ;
        if (++cnt >= 20){
            cout << "..." << endl;
            break;
        }
    }
    cout << endl;
};

double DecisionTree::findCutoffValue(string gene, list<unsigned int> sampleList){
    /* Goal of this method (function):
    // Finding a cutoff value should be implemented and the cutoff value should be returned.
    // 1. For samples in `sampleList`, get expression value of a gene.
    // 2. Within the samples in `sampleList`, find the cutoff value of expression of a gene.
    // 3. return the cutoff value (double)
    */

    list<unsigned int>::iterator it;
    double SuM_norm = 0; double SuM_cancr = 0; int cnt_n = 0; int cnt_c = 0;
    for(it = sampleList.begin(); it!=sampleList.end(); it++){
        if (sampleCls.find(*it)->second) {
            SuM_norm += exp.find(pair<unsigned int, string>(*it, gene))->second;
            cnt_n++;
        }
        else if (!(sampleCls.find(*it)->second)) {
            SuM_cancr += exp.find(pair<unsigned int, string>(*it, gene))->second;
            cnt_c++;
        }
    }
    if ( cnt_n == 0){
        return (SuM_cancr / (double)cnt_c);
    }else if ( cnt_c == 0 ){
        return (SuM_norm / (double)cnt_n);
    }

    return ((SuM_norm / (double)cnt_n) + (SuM_cancr / (double)cnt_c)) / 2;
};
void DecisionTree::testCutoff(string gene, list<unsigned int> sampleList){
    /* Test4
	// This is for testing your code.
	// Printing out the classes and indeces of all samples in `sampleList` */
    int cnt = 0;
    cout << "TEST4: print class of samples in 'sampleList'" << endl;
    list<unsigned int>::iterator it;
    for(it= sampleList.begin(); it!=sampleList.end(); it++){
        cout << "sample#: " << *it << "\t";
        cout << "sample# in map: " << sampleCls.find(*it)->first << "\t"; // get sample number
        cout << "sample class in map: " << (sampleCls.find(*it)->second  ? "Normal" : "Disease") << endl; // get corresponding sample class
        if (++cnt >= 20){
            cout << "..." << endl;
            break;
        }
    }
    cout << endl;

    /* Test5
    // This is for testing your code.
    // Printing out the expression value of a gene for all samples in `sampleList` */
    cnt = 0;
    cout << "TEST5: find expression value for gene of given sample" << endl;
    cout << "input gene: " << gene << endl;
    for(it= sampleList.begin(); it!=sampleList.end(); it++){
        cout << "Expression value of " << gene << " for sample " << *it << ": "; // print gene names
        cout << exp.find(pair<unsigned int, string>(*it, gene))->second << endl; // find value(expression value) for sample number and gene
        if (++cnt >= 20){
            cout << "..." << endl;
            break;
        }
    }
    cout << endl;

    /* Test6
    // This is for test your method(function `findCutoffValue`).
    // Print out the cutoff value for a certain gene
    */
    cout << "TEST6: Find cutoff value for a gene of given samples" << endl;
    cout << "input gene: " << gene << endl;
    cout << "cutoffvalue : " << findCutoffValue(gene, sampleList) << endl;
    cout << endl;
};

double DecisionTree::getInfoGain(string gene, double cutoffValue, list<unsigned int> sampleList){
    /* Goal of this method (function):
    // For given samples and gene, calculating information gain part should be implemented and the value is returned.
    // 1. For samples in `sampleList`, get expression value of a gene.
    // 2. Within the samples in `sampleList`, calculate information gain of separation of sample set with a certain cutoff value.
    // 3. return the cutoff value (double)
    */

    // TIP! : Information gain should be calculated ONLY for inputted samples.

    double entro_1, entro_h, entro_l;
    int Norm = 0; int Cancr = 0; int Norm_h = 0; int Norm_l = 0; int Cancr_h = 0; int Cancr_l = 0;
    list<unsigned int>::iterator it;
    for (it = sampleList.begin(); it != sampleList.end(); it++) {
        if (sampleCls.find(*it)->second) {
            Norm++;
            if (exp.find(pair<unsigned int, string>(*it, gene))->second >= cutoffValue) Norm_h++;
            else Norm_l++;
        }
        else if (!(sampleCls.find(*it)->second)){
            Cancr++;
            if (exp.find(pair<unsigned int, string>(*it, gene))->second >= cutoffValue) Cancr_h++;
            else Cancr_l++;
        }
    }
    entro_1 = entropycal(Norm, Cancr);
    entro_h = entropycal(Norm_h, Cancr_h);
    entro_l = entropycal(Norm_l, Cancr_l);

    double Output = entro_1 - (((double)Norm_h + (double)Cancr_h) / ((double)Norm + (double)Cancr))*entro_h
                    - (((double)Norm_l + (double)Cancr_l) / ((double)Norm + (double)Cancr))*entro_l;

    return Output;
};
void DecisionTree::testInfoGain(string gene, list<double> cutoffValueList, list<unsigned int> sampleList){
    /* Test6
	// This is for testing your method(function `getInfoGain`).
	// Printing out the InfoGain value for expression value of a gene for all samples in `sampleList` */
    list<double>::iterator it;
    cout << "TEST7: Check information gain value for a gene of given samples with various cutoff value" << endl;
    cout << "input gene: " << gene << endl;
    for(it= cutoffValueList.begin(); it!=cutoffValueList.end(); it++){
        cout << "InfoGain of " << gene << " for samples with cutoff " << *it <<" : "; // print gene name
        cout << getInfoGain(gene, *it, sampleList) << endl; // find value(expression value) for sample number and gene
    }
    cout << endl;
}


int main()
{
    DecisionTree d;
    d.readDataFromFile("data/madata.txt"); // change to proper filename for your case
    list<string> deg = d.readDEGFromFile("data/deg.txt"); // change to proper filename for your case
    d.testDEGList(deg);

    list<unsigned int> l;
    l.push_back(1);
    l.push_back(2);
    l.push_back(4);
    l.push_back(5);
    l.push_back(6);
    l.push_back(3);

    list<string>::iterator it = deg.begin();
    double tmp_cov = d.findCutoffValue(*it,l);
    d.testCutoff(*it, l);
    cout << tmp_cov << endl;

    list<double> cov_l;
    cov_l.push_back(tmp_cov);
    cov_l.push_back(3.1);
    cov_l.push_back(4.4);
    d.testInfoGain(*it, cov_l, l);
    // system("PAUSE");
}
