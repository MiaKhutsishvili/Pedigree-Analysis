/*
    Pedigree Analysis
    Author: Kimia Khoodsiyani
*/

#include <iostream>
#include <fstream>
#include <vector>
#include <queue>
#include <cmath>
#include <algorithm>
#include <string>

using namespace std;

// People
const int maxNumOfPeople = 2e5;
int numOfPeople;
vector<int> gender; // 1 if Female.
vector<bool> AffectionFile;

// Pedigree
int Partner[maxNumOfPeople];
vector<pair<int, int>> Parents;
vector<int> Children[maxNumOfPeople];

// bfs
queue<int> Headers;
bool bfsMarker[maxNumOfPeople];
bool subTreeMarker[maxNumOfPeople]; // Flaged with the mother

// Probability
double Gene_Prob[maxNumOfPeople][5][7];
double Model_Prob[5];
/*
    5 Modes
    Female Genes: {DD, DR, RR, XdXd, XdXr, XrXr}
    Male Genes: {DD, DR, RR, XdY, XrY, XY*, XY}
*/
pair<double, string> mostLikableModel[5];

void bfs(int);
void autoDomProb(int);
void autoRecProb(int);
void xLinkedDomProb(int);
void xLinkedRecProb(int);
void yLinked(int);
void PedStarters(int);

double ChildrenAffectionProb(int, double);
double DaughtersAffectionProb(int, double);
double SonsAffectionProb(int, double);
double ADModelProb(int);
double ARModelProb(int);
double XLDModeProbe(int);
double XLRModeProbe(int);
double YLModeProbe(int);


void exportToDOT(const string &filename)
{
    ofstream out(filename);
    out << "digraph Pedigree {\n";
    out << "  rankdir=TB;\n"; // Tree layout: top to bottom
    out << "  node [fontname=\"Arial\", fontsize=12, style=filled, fillcolor=white];\n";
    out << "  edge [color=gray50];\n\n";

    // Draw individuals
    for (int i = 0; i < numOfPeople; i++)
    {
        string nodeName = "Person_" + to_string(i + 1);
        string shape = gender[i] ? "circle" : "box"; // Female: circle, Male: box
        string color = AffectionFile[i] ? "red" : "black";
        string genderSymbol = gender[i] ? "♀" : "♂";

        out << "  " << nodeName << " [label=\"" << genderSymbol << " " << i + 1
            << "\", shape=" << shape << ", color=" << color << "];\n";
    }

    out << "\n";

    // Draw couples and children
    for (int i = 0; i < numOfPeople; i++)
    {
        int partner = Partner[i];
        if (partner != -1 && i < partner)
        {
            string coupleNode = "Couple_" + to_string(i + 1) + "_" + to_string(partner + 1);

            // Align couple horizontally
            out << "  { rank=same; Person_" << i + 1 << "; Person_" << partner + 1 << "; }\n";

            // Invisible anchor node
            out << "  " << coupleNode << " [shape=point, width=0, label=\"\"];\n";

            // Horizontal line between partners
            out << "  Person_" << i + 1 << " -> " << coupleNode << " [dir=none, constraint=false, color=gray70];\n";
            out << "  Person_" << partner + 1 << " -> " << coupleNode << " [dir=none, constraint=false, color=gray70];\n";

            // Children descend from couple node
            for (int child : Children[i])
            {
                if (Parents[child].first == i || Parents[child].second == i)
                    out << "  " << coupleNode << " -> Person_" << child + 1 << " [color=gray30];\n";
            }

            out << "\n";
        }
    }

    out << "}\n";
    out.close();
}

int main()
{
    for (int i = 0; i < 5; i++)
        Model_Prob[i] = 1.0;

    // Getting the input pedigree...
    cout << "Hello dear user!" << endl;
    cout << "Please enter the number of people in your pedigree: ";
    cin >> numOfPeople;
    for (int i = 0; i < numOfPeople; i++)
    {
        bool gndr;
        cout << "Enter the gender of person #" << i + 1 << ":" << endl;
        cout << "1 if Female / 0 if Male -> ";
        cin >> gndr;
        gender.push_back(gndr);

        pair<int, int> prnts;
        cout << "Enter parents of #" << i + 1 << ":" << endl;
        cout << "No parents? Just enter 0 as their parents." << endl;
        cout << "#" << i + 1 << "'s mother -> ";
        cin >> prnts.first;
        cout << "#" << i + 1 << "'s father -> ";
        cin >> prnts.second;
        prnts.first--;
        prnts.second--;
        Parents.push_back(prnts);
        if (prnts.first != -1)
        {
            Children[prnts.first].push_back(i);
            Children[prnts.second].push_back(i);
        }
        else
            Headers.push(i);

        cout << "Enter the partner of #" << i + 1 << ":" << endl;
        cout << "No partner? Enter 0." << endl;
        cout << "-> ";
        cin >> Partner[i];
        Partner[i]--;

        bool Affected;
        cout << "Is #" << i + 1 << " affected? " << endl;
        cout << "1 if yes / 0 if not -> ";
        cin >> Affected;
        AffectionFile.push_back(Affected);
    }

    exportToDOT("pedigree.dot");

    // Gene Probability Calculation
    while (Headers.size())
        bfs(Headers.front());

    // Best Model
    /* Normilizing:
    double total = 0;
    for (int i = 0; i < 5; i++)
        total += Model_Prob[i]; // Assuming uniform priors

    if (total > 0)
    {
        for (int i = 0; i < 5; i++)
            Model_Prob[i] = Model_Prob[i] / total;
    }
    */
    mostLikableModel[0].second = "Autosome Dominant Inheritance";
    mostLikableModel[1].second = "Autosome Recessive Inheritance";
    mostLikableModel[2].second = "X_Linked Dominant Inheritance";
    mostLikableModel[3].second = "X_Linked Recessive Inheritance";
    mostLikableModel[4].second = "Y_Linked Inheritance";
    for (int i = 0; i < 5; i++)
        mostLikableModel[i].first = Model_Prob[i];
    sort(mostLikableModel, mostLikableModel + 5);
    reverse(mostLikableModel, mostLikableModel + 5);
    cout << "Best Model is : " << mostLikableModel[0].second << endl;
    cout << "The lilkelihood of each model is:" << endl;
    for (int i = 0; i < 5; i++)
    {
        cout << i + 1 << ") ";
        cout << mostLikableModel[i].first << " for " << mostLikableModel[i].second << endl;
    }

    return 0;
}

void bfs(int node)
{
    if (Parents[node].first == -1)
    {
        PedStarters(node);
        bfsMarker[node] = 1;
    }

    int spouse = Partner[node];
    int mother;
    if (gender[node])
        mother = node;
    else
        mother = spouse;

    if (spouse != -1 && bfsMarker[spouse] && !subTreeMarker[mother])
    {
        subTreeMarker[mother] = 1;
        Model_Prob[0] *= ADModelProb(mother);
        Model_Prob[1] *= ARModelProb(mother);
        Model_Prob[2] *= XLDModeProbe(mother);
        Model_Prob[3] *= XLRModeProbe(mother);
        Model_Prob[4] *= YLModeProbe(mother);

        for (int i = 0; i < Children[node].size(); i++)
        {
            int Child = Children[node][i];
            autoDomProb(Child);
            autoRecProb(Child);
            xLinkedDomProb(Child);
            xLinkedRecProb(Child);
            yLinked(Child);
            bfsMarker[Child] = 1;
            Headers.push(Child);
        }
    }
    // Else? We will get here from the other parent!

    Headers.pop();
    return;
}

void PedStarters(int node)
{
    if (!Children[node].size())
    {
        Model_Prob[0] *= ((AffectionFile[node] * 2.0 / 3.0) + (!AffectionFile[node] * 1.0 / 3.0));
        Model_Prob[1] *= ((AffectionFile[node] * 1.0 / 3.0) + (!AffectionFile[node] * 2.0 / 3.0));
        if (gender[node])
        {
            Model_Prob[2] *= ((AffectionFile[node] * 2.0 / 3.0) + (!AffectionFile[node] * 1.0 / 3.0));
            Model_Prob[3] *= ((AffectionFile[node] * 1.0 / 3.0) + (!AffectionFile[node] * 2.0 / 3.0));
            Model_Prob[4] *= ((AffectionFile[node] * 0) + (!AffectionFile[node] * 1));
        }
        else
        {
            Model_Prob[2] *= ((AffectionFile[node] * 1.0 / 2.0) + (!AffectionFile[node] * 1.0 / 2.0));
            Model_Prob[3] *= ((AffectionFile[node] * 1.0 / 2.0) + (!AffectionFile[node] * 1.0 / 2.0));
            Model_Prob[4] *= ((AffectionFile[node] * 1.0 / 2.0) + (!AffectionFile[node] * 1.0 / 2.0));
        }
    }
    if (AffectionFile[node])
    {
        // Mode 0 -> DD or DR
        Gene_Prob[node][0][0] = 0.5;
        Gene_Prob[node][0][1] = 0.5;

        // Mode 1 -> RR
        Gene_Prob[node][1][2] = 1;

        // Mode 2
        if (gender[node])
        {
            // XdXd or XdXr
            Gene_Prob[node][2][3] = 0.5;
            Gene_Prob[node][2][4] = 0.5;
        }
        else
            // XdY
            Gene_Prob[node][2][3] = 1;

        // Mode 3
        if (gender[node])
            // XrXr
            Gene_Prob[node][3][5] = 1;
        else
            // XrY
            Gene_Prob[node][3][4] = 1;

        // Mode 4
        if (!gender[node])
            // XY*
            Gene_Prob[node][4][5] = 1;

        // Others are 0 by default.
    }
    else
    {
        // Mode 0 -> RR
        Gene_Prob[node][0][2] = 1;

        // Mode 1 -> DD or DR
        Gene_Prob[node][1][0] = 0.5;
        Gene_Prob[node][1][1] = 0.5;

        // Mode 2
        if (gender[node])
            // XrXr
            Gene_Prob[node][2][5] = 1;
        else
            // XrY
            Gene_Prob[node][2][4] = 1;

        // Mode 3
        if (gender[node])
        {
            // XdXd or XdXr
            Gene_Prob[node][3][3] = 0.5;
            Gene_Prob[node][3][4] = 0.5;
        }
        else
            // XdY
            Gene_Prob[node][3][3] = 1;

        // Mode 4
        if (!gender[node])
            // XY
            Gene_Prob[node][4][6] = 1;

        // Others are 0 by default.
    }

    return;
}

void autoDomProb(int node)
{
    int mom = Parents[node].first;
    int dad = Parents[node].second;

    // DD / Affected.
    Gene_Prob[node][0][0] = (Gene_Prob[mom][0][0] + 0.5 * Gene_Prob[mom][0][1]) *
                            (Gene_Prob[dad][0][0] + 0.5 * Gene_Prob[dad][0][1]);

    // RR
    Gene_Prob[node][0][2] = (Gene_Prob[mom][0][2] + 0.5 * Gene_Prob[mom][0][1]) *
                            (Gene_Prob[dad][0][2] + 0.5 * Gene_Prob[dad][0][1]);

    // DR / Affected.
    Gene_Prob[node][0][1] = 1 - (Gene_Prob[node][0][0] + Gene_Prob[node][0][2]);

    // Considering the facts...
    Gene_Prob[node][0][0] *= AffectionFile[node];
    Gene_Prob[node][0][1] *= AffectionFile[node];
    Gene_Prob[node][0][2] *= !AffectionFile[node];

    return;
}

void autoRecProb(int node)
{
    int mom = Parents[node].first;
    int dad = Parents[node].second;

    // DD
    Gene_Prob[node][1][0] = (Gene_Prob[mom][1][0] + 0.5 * Gene_Prob[mom][1][1]) *
                            (Gene_Prob[dad][1][0] + 0.5 * Gene_Prob[dad][1][1]);

    // RR Affected.
    Gene_Prob[node][1][2] = (Gene_Prob[mom][1][2] + 0.5 * Gene_Prob[mom][1][1]) *
                            (Gene_Prob[dad][1][2] + 0.5 * Gene_Prob[dad][1][1]);

    // DR
    Gene_Prob[node][1][1] = 1 - (Gene_Prob[node][1][0] + Gene_Prob[node][1][2]);

    // Considering the facts...
    Gene_Prob[node][1][0] *= !AffectionFile[node];
    Gene_Prob[node][1][1] *= !AffectionFile[node];
    Gene_Prob[node][1][2] *= AffectionFile[node];

    return;
}

void xLinkedDomProb(int node)
{
    int mom = Parents[node].first;
    int dad = Parents[node].second;

    if (gender[node])
    {
        // XdXd Affected.
        Gene_Prob[node][2][3] = Gene_Prob[dad][2][3] *
                                (Gene_Prob[mom][2][3] + 0.5 * Gene_Prob[mom][2][4]);

        // XrXr
        Gene_Prob[node][2][5] = Gene_Prob[dad][2][4] *
                                (Gene_Prob[mom][2][5] + 0.5 * Gene_Prob[mom][2][4]);

        // XrXd Affected.
        Gene_Prob[node][2][4] = 1 - (Gene_Prob[node][2][3] + Gene_Prob[node][2][5]);

        // Considering the facts...
        Gene_Prob[node][2][3] *= AffectionFile[node];
        Gene_Prob[node][2][4] *= AffectionFile[node];
        Gene_Prob[node][2][5] *= !AffectionFile[node];
    }
    else
    {
        // XdY Affected.
        Gene_Prob[node][2][3] = Gene_Prob[mom][2][3] + 0.5 * Gene_Prob[mom][2][4];

        // XrY
        Gene_Prob[node][2][4] = Gene_Prob[mom][2][5] + 0.5 * Gene_Prob[mom][2][4];

        // Considering the facts...
        Gene_Prob[node][2][3] *= AffectionFile[node];
        Gene_Prob[node][2][4] *= !AffectionFile[node];
    }

    return;
}

void xLinkedRecProb(int node)
{
    int mom = Parents[node].first;
    int dad = Parents[node].second;

    if (gender[node])
    {
        // XdXd
        Gene_Prob[node][3][3] = Gene_Prob[dad][3][3] *
                                (Gene_Prob[mom][3][3] + 0.5 * Gene_Prob[mom][3][4]);

        // XrXr Affected.
        Gene_Prob[node][3][5] = Gene_Prob[dad][3][4] *
                                (Gene_Prob[mom][3][5] + 0.5 * Gene_Prob[mom][3][4]);

        // XrXd
        Gene_Prob[node][3][4] = 1 - (Gene_Prob[node][3][3] + Gene_Prob[node][3][5]);

        // Considering the facts...
        Gene_Prob[node][3][3] *= !AffectionFile[node];
        Gene_Prob[node][3][4] *= !AffectionFile[node];
        Gene_Prob[node][3][5] *= AffectionFile[node];
    }
    else
    {
        // XdY
        Gene_Prob[node][3][3] = Gene_Prob[mom][3][3] + 0.5 * Gene_Prob[mom][3][4];

        // XrY Affected.
        Gene_Prob[node][3][4] = Gene_Prob[mom][3][5] + 0.5 * Gene_Prob[mom][3][4];

        // Considering the facts...
        Gene_Prob[node][3][3] *= !AffectionFile[node];
        Gene_Prob[node][3][4] *= AffectionFile[node];
    }

    return;
}

void yLinked(int node)
{
    int mom = Parents[node].first;
    int dad = Parents[node].second;

    if (gender[node])
        return;
    else
    {
        // XY* Affected.
        Gene_Prob[node][4][5] = Gene_Prob[dad][4][5];

        // XY
        Gene_Prob[node][4][6] = Gene_Prob[dad][4][6];

        // Considering the facts...
        Gene_Prob[node][4][5] *= AffectionFile[node];
        Gene_Prob[node][4][6] *= !AffectionFile[node];
    }

    return;
}

double ChildrenAffectionProb(int mother, double AffProb)
{
    double prob = 1.0;
    for (int i = 0; i < Children[mother].size(); i++)
    {
        int Child = Children[mother][i];
        prob *= ((AffectionFile[Child] * AffProb) + (!AffectionFile[Child] * (1 - AffProb)));
    }
    return prob;
}

double DaughtersAffectionProb(int mother, double AffProb)
{
    double prob = 1.0;
    for (int i = 0; i < Children[mother].size(); i++)
    {
        int Child = Children[mother][i];
        if (gender[Child])
            prob *= ((AffectionFile[Child] * AffProb) + (!AffectionFile[Child] * (1 - AffProb)));
    }
    return prob;
}

double SonsAffectionProb(int mother, double AffProb)
{
    double prob = 1.0;
    for (int i = 0; i < Children[mother].size(); i++)
    {
        int Child = Children[mother][i];
        if (!gender[Child])
            prob *= ((AffectionFile[Child] * AffProb) + (!AffectionFile[Child] * (1 - AffProb)));
    }
    return prob;
}

double ADModelProb(int mother)
{
    int father = Partner[mother];
    double Possibility[9];
    double answer = 0;

    // DD , DD -> All Children Affected.
    Possibility[0] = Gene_Prob[mother][0][0] * Gene_Prob[father][0][0];
    Possibility[0] *= ChildrenAffectionProb(mother, 1);

    // DD , DR -> All Children Affected.
    Possibility[1] = Gene_Prob[mother][0][0] * Gene_Prob[father][0][1];
    Possibility[1] *= ChildrenAffectionProb(mother, 1);

    // DD , RR -> All Children Affected.
    Possibility[2] = Gene_Prob[mother][0][0] * Gene_Prob[father][0][2];
    Possibility[2] *= ChildrenAffectionProb(mother, 1);

    // DR , DD -> All Children Affected.
    Possibility[3] = Gene_Prob[mother][0][1] * Gene_Prob[father][0][0];
    Possibility[3] *= ChildrenAffectionProb(mother, 1);

    // DR , DR -> 3/4 Of Children Affected.
    Possibility[4] = Gene_Prob[mother][0][1] * Gene_Prob[father][0][1];
    Possibility[4] *= ChildrenAffectionProb(mother, 0.75);

    // DR , RR -> 1/2 Of Chidren Affected.
    Possibility[5] = Gene_Prob[mother][0][1] * Gene_Prob[father][0][2];
    Possibility[5] *= ChildrenAffectionProb(mother, 0.5);

    // RR , DD -> All Chidren Affected.
    Possibility[6] = Gene_Prob[mother][0][2] * Gene_Prob[father][0][0];
    Possibility[6] *= ChildrenAffectionProb(mother, 1);

    // RR , DR -> 1/2 Of Chidren Affected.
    Possibility[7] = Gene_Prob[mother][0][2] * Gene_Prob[father][0][1];
    Possibility[7] *= ChildrenAffectionProb(mother, 0.5);

    // RR , RR -> None Of Chidren Affected.
    Possibility[8] = Gene_Prob[mother][0][2] * Gene_Prob[father][0][2];
    Possibility[8] *= ChildrenAffectionProb(mother, 0);

    for (int i = 0; i < 9; i++)
        answer += Possibility[i];
    return answer;
}

double ARModelProb(int mother)
{
    int father = Partner[mother];
    double Possibility[9];
    double answer = 0;

    // DD , DD -> None of Children Affected.
    Possibility[0] = Gene_Prob[mother][1][0] * Gene_Prob[father][1][0];
    Possibility[0] *= ChildrenAffectionProb(mother, 0);

    // DD , DR -> None of Children Affected.
    Possibility[1] = Gene_Prob[mother][1][0] * Gene_Prob[father][1][1];
    Possibility[1] *= ChildrenAffectionProb(mother, 0);

    // DD , RR -> None of Children Affected.
    Possibility[2] = Gene_Prob[mother][1][0] * Gene_Prob[father][1][2];
    Possibility[2] *= ChildrenAffectionProb(mother, 0);

    // DR , DD -> None of Children Affected.
    Possibility[3] = Gene_Prob[mother][1][1] * Gene_Prob[father][1][0];
    Possibility[3] *= ChildrenAffectionProb(mother, 0);

    // DR , DR -> 1/4 Of Children Affected.
    Possibility[4] = Gene_Prob[mother][1][1] * Gene_Prob[father][1][1];
    Possibility[4] *= ChildrenAffectionProb(mother, 0.25);

    // DR , RR -> 1/2 Of Chidren Affected.
    Possibility[5] = Gene_Prob[mother][1][1] * Gene_Prob[father][1][2];
    Possibility[5] *= ChildrenAffectionProb(mother, 0.5);

    // RR , DD -> None of Chidren Affected.
    Possibility[6] = Gene_Prob[mother][1][2] * Gene_Prob[father][1][0];
    Possibility[6] *= ChildrenAffectionProb(mother, 0);

    // RR , DR -> 1/2 Of Chidren Affected.
    Possibility[7] = Gene_Prob[mother][1][2] * Gene_Prob[father][1][1];
    Possibility[7] *= ChildrenAffectionProb(mother, 0.5);

    // RR , RR -> All Chidren Affected.
    Possibility[8] = Gene_Prob[mother][1][2] * Gene_Prob[father][1][2];
    Possibility[8] *= ChildrenAffectionProb(mother, 1);

    for (int i = 0; i < 9; i++)
        answer += Possibility[i];
    return answer;
}

double XLDModeProbe(int mother)
{
    int father = Partner[mother];
    double Possibility[6];
    double answer = 0;

    // XdXd , XdY -> ALl Children Affected.
    Possibility[0] = Gene_Prob[mother][2][3] * Gene_Prob[father][2][3];
    Possibility[0] *= ChildrenAffectionProb(mother, 1);

    // XdXd , XrY -> All Children Affected.
    Possibility[1] = Gene_Prob[mother][2][3] * Gene_Prob[father][2][4];
    Possibility[1] *= ChildrenAffectionProb(mother, 1);

    // XdXr , XdY -> All Daughters Affected,
    //               1/2 of Sons Affected.
    Possibility[2] = Gene_Prob[mother][2][4] * Gene_Prob[father][2][3];
    Possibility[2] *= DaughtersAffectionProb(mother, 1);
    Possibility[2] *= SonsAffectionProb(mother, 0.5);

    // XdXr , XrY -> 1/2 of Daughters Affected,
    //               1/2 of Sons Affected.
    Possibility[3] = Gene_Prob[mother][2][4] * Gene_Prob[father][2][4];
    Possibility[3] *= DaughtersAffectionProb(mother, 0.5);
    Possibility[3] *= SonsAffectionProb(mother, 0.5);
    // XrXr , XdY -> All Daughters Affected,
    //               None of Sons Affected.
    Possibility[4] = Gene_Prob[mother][2][5] * Gene_Prob[father][2][3];
    Possibility[4] *= DaughtersAffectionProb(mother, 1);
    Possibility[4] *= SonsAffectionProb(mother, 0);

    // XrXr , XrY -> None of Daughters Affected,
    //               None of Sons Affected.
    Possibility[5] = Gene_Prob[mother][2][5] * Gene_Prob[father][2][4];
    Possibility[5] *= DaughtersAffectionProb(mother, 0);
    Possibility[5] *= SonsAffectionProb(mother, 0);

    for (int i = 0; i < 6; i++)
        answer += Possibility[i];
    return answer;
}

double XLRModeProbe(int mother)
{
    int father = Partner[mother];
    double Possibility[6];
    double answer = 0;

    // XdXd , XdY -> None of Children Affected.
    Possibility[0] = Gene_Prob[mother][3][3] * Gene_Prob[father][3][3];
    Possibility[0] *= ChildrenAffectionProb(mother, 0);

    // XdXd , XrY -> None of Children Affected.
    Possibility[1] = Gene_Prob[mother][3][3] * Gene_Prob[father][3][4];
    Possibility[1] *= ChildrenAffectionProb(mother, 0);

    // XdXr , XdY -> None of Daughters Affected,
    //               1/2 of Sons Affected.
    Possibility[2] = Gene_Prob[mother][3][4] * Gene_Prob[father][3][3];
    Possibility[2] *= DaughtersAffectionProb(mother, 0);
    Possibility[2] *= SonsAffectionProb(mother, 0.5);

    // XdXr , XrY -> 1/2 of Daughters Affected,
    //               1/2 of Sons Affected.
    Possibility[3] = Gene_Prob[mother][3][4] * Gene_Prob[father][3][4];
    Possibility[3] *= DaughtersAffectionProb(mother, 0.5);
    Possibility[3] *= SonsAffectionProb(mother, 0.5);

    // XrXr , XdY -> None of Daughters Affected,
    //               ALl Sons Affected.
    Possibility[4] = Gene_Prob[mother][3][5] * Gene_Prob[father][3][3];
    Possibility[4] *= DaughtersAffectionProb(mother, 0);
    Possibility[4] *= SonsAffectionProb(mother, 1);

    // XrXr , XrY -> All Children Affected.
    Possibility[5] = Gene_Prob[mother][3][5] * Gene_Prob[father][3][4];
    Possibility[5] *= ChildrenAffectionProb(mother, 1);

    for (int i = 0; i < 6; i++)
        answer += Possibility[i];
    return answer;
}

double YLModeProbe(int mother)
{
    int father = Partner[mother];
    double Possibility[2];

    // No Woman Has it.
    if (AffectionFile[mother])
        return 0;
    double answer = 0;

    // XY* -> All Sons Affected.
    Possibility[0] = Gene_Prob[father][4][5];
    Possibility[0] *= SonsAffectionProb(mother, 1);

    // XY -> None of Sons Affected.
    Possibility[1] = Gene_Prob[father][4][6];
    Possibility[1] *= SonsAffectionProb(mother, 0);

    for (int i = 0; i < 2; i++)
        answer += Possibility[i];

    answer *= DaughtersAffectionProb(mother, 0);

    return answer;
}