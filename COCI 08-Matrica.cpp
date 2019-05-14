/*
Notice that if we make the entire matrix it will be of [30000][30000] size. There is a reason why <= 50 columns have to be output

The answer is impossible if number of odd frequency characters is greater than diagonal length i.e. n (Why?) Because all odd frequency elements have to have at least 1 element in the diagonal else there will be at least 1 character whose pair doesn't exist.

Let's say we are solving for row number 'row'. Notice that the matrix is symmetric. So if we put X at matrix[row][col] then matrix[col][row] has to be X as well
So we just solve for the upper triangular matrix plus the diagonal. So at row 'row' we start from cell matrix[row][row].

We break solving each row into 2 parts.

1. Solve the diagonal element
	We count the number of odd frequency characters (we have to do it in every iteration since frequency changes)
	if N - row + 1 (length of diagonal from [row][row] to [n][n]) is greater than the number of odd elements we can afford to place a lexicographically smaller element (may be even or odd doesnt matter)
	but if it is equal than we have to place the smallest odd frequency element.

2. Solve for the remaining of the row
	If the diagonal element [row][row] was in a required column, then we have to fill that column completely (mat[][row]), so we iterate over all the elements on the row and set the smallest character if it is >= 2 since we have to place it in the symmetric position as well
	If the diagonal is not a required column, we try to iterate over the required columns in that row to save time. We know how many same characters we can place is freq[character]/2 so we directly go to next required column or till we have to change the character because the frequency of that character became 0

Note that we solve the matrix row by row so after each iteration on row we have the answer for that row.	
*/

#include <bits/stdc++.h>
using namespace std;
#define endl '\n'

char ch;
char matrix[30005][55];
int x, n, k, idx, odd, p[55], P, freq[26];

int main(){
    ios_base::sync_with_stdio(0);cin.tie(NULL);cout.tie(NULL);
	
    cin >> n >> k;

    for(int i = 0; i < k; ++i){
    	cin >> ch >> x;
    	freq[ch - 'A'] = x;
    	if(x & 1)odd++;
    }

    cin >> P;
    for(int i = 0; i < P; ++i)
    	cin >> p[i];
    p[P] = n + 1;

    if(odd > n){
    	cout << "IMPOSSIBLE" << endl;
    	return 0;
    }

    int p_start = 0, character = 0;

    for(int row = 1; row <= n; ++row){

    	int first = 26, first_odd = 26;

    	odd = 0;
    	for(int i = 25; i >= 0; --i){
    		if(freq[i] > 0)first = i;
    		if(freq[i] & 1){first_odd = i; ++odd;}
    	}

    	if(n - row + 1 > odd)character = first; //can place 2 even ones (if they are smaller) because number of odd ones are less than diagonal
    	else character = first_odd;

    	//handle diagonal
    	freq[character]--; //place on diagonal
    	int required_column = (row == p[p_start]);
    	if(required_column)matrix[row][p_start] = char(character + 'A');

    	//handle rest of row
    	int iterator = row + 1;
    	character = 0;
    	int position = p_start + required_column;

    	while(iterator <= n){
    		while(freq[character] < 2)character++; //have to have atleast 2 to place
    		if(required_column)matrix[iterator][p_start] = char(character + 'A');
    		if(iterator == p[position])matrix[row][position++] = char(character + 'A');
    		int next_iterator = p[position];
    		if(required_column)next_iterator = iterator + 1;
    		next_iterator = min(next_iterator, iterator + freq[character]/2);
    		freq[character] = freq[character] - (next_iterator - iterator)*2;
    		iterator = next_iterator;
    	}

    	if(required_column)
    		p_start++;

    	for(int col = 0; col < P; ++col)
    		cout << matrix[row][col];
    	cout << endl;

    }

    return 0;
}

