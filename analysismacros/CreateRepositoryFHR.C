#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include <vector>
#include <TCanvas.h>
#include <TGraph.h>
#include <TFile.h>
#include <TH1.h>
#include <TROOT.h>
#include <TSystem.h>


using namespace std;
//function to create the right layout for comments in the terminal output
void extended_to_reduced_terminal(vector<Int_t> extended, vector<Int_t> reduced)
{
	Int_t dimension_extended, dimension_reduced;
	Int_t low_edge_reduced, high_edge_reduced, width_reduced;
	Int_t low_edge_extended = 0;
	Int_t high_edge_extended = 0;
	Int_t width_extended = 0;
	dimension_extended = extended.size();
	reduced.push_back(extended[0]);
	for(Int_t i=1; i < dimension_extended; i++)
	{
		if (extended[i]!=extended[i-1]+1 || extended[i]!=extended[i+1]-1) reduced.push_back(extended[i]);
	}
						
	dimension_reduced = reduced.size();
			
	for(Int_t i=0; i < dimension_reduced-1; i++)
	{
		if(reduced[i+1]-reduced[i]==1) {cout << reduced[i] << ", ";}
		else 
		{
			low_edge_reduced = reduced[i];
			high_edge_reduced = reduced[i+1];
			width_reduced = high_edge_reduced - low_edge_reduced; 
			for(Int_t i=0; i < dimension_extended; i++)
			{
				if(extended[i]==low_edge_reduced) low_edge_extended = i;
				if(extended[i]==high_edge_reduced) high_edge_extended = i;
			}
			width_extended = high_edge_extended -low_edge_extended;
			if(width_reduced==width_extended) {cout << reduced[i] << "-";} 
			else {cout << reduced[i] << ", "; } 
		}
	}
	cout << reduced[dimension_reduced-1];
};

//function to create the right layout for comments in the repository
void extended_to_reduced_repository(vector<Int_t> extended, vector<Int_t> reduced, string& build_comment)
{
	Int_t dimension_extended, dimension_reduced;
	Int_t low_edge_reduced, high_edge_reduced, width_reduced;
	Int_t low_edge_extended = 0; 
	Int_t high_edge_extended = 0;
	Int_t width_extended = 0;
	dimension_extended = extended.size();
	reduced.push_back(extended[0]);
	for(Int_t i=1; i < dimension_extended; i++)
	{
		if (extended[i]!=extended[i-1]+1 || extended[i]!=extended[i+1]-1) reduced.push_back(extended[i]);
	}
						
	dimension_reduced = reduced.size();
	string convert = "";
			
	for(Int_t i=0; i < dimension_reduced-1; i++)
	{
		if(reduced[i+1]-reduced[i]==1) 
		{
			convert = to_string(reduced[i]);
			build_comment.append(convert);
			build_comment.append(", ");
		}
		else 
		{
			low_edge_reduced = reduced[i];
			high_edge_reduced = reduced[i+1];
			width_reduced = high_edge_reduced - low_edge_reduced; 
			for(Int_t i=0; i < dimension_extended; i++)
			{
				if(extended[i]==low_edge_reduced) low_edge_extended = i;
				if(extended[i]==high_edge_reduced) high_edge_extended = i;
			}
			width_extended = high_edge_extended -low_edge_extended;
			if(width_reduced==width_extended) 
			{
				convert = to_string(reduced[i]);
				build_comment.append(convert);
				build_comment.append("-");
			} 
			else 
			{
				convert = to_string(reduced[i]);
				build_comment.append(convert);
				build_comment.append(", ");
			} 
		}
	}
	convert = to_string(reduced[dimension_reduced-1]);
	build_comment.append(convert);
};

//analysis for the terminal output
void terminal_output(Int_t layer, Int_t n_stave, Int_t n_events, vector<string> cumulative_label, vector<string> label_event, vector<TGraph*> staves, const char* layer_name)
{
	//definition of counters to help the division
	Int_t counter_inactive = 0;
	Int_t counter_less8 = 0; //counter for FHR < 10^-8
	Int_t counter_less4 = 0; //counter for 10^-8 < FHR < 10^-4 
	Int_t counter_less3 = 0; //counter for 10^-4 < FHR < 10^-3
	Int_t counter_more3 = 0; //counter for FHR > 10^-3
	Int_t counter_emulated_pp = 0; //counter for emulated pp
	Int_t counter_emulated_PbPb = 0; //counter for emulated PbPb
	
	//definition of vector<int> to store the different behaviour extended
	vector<Int_t> staves_inactive;
	vector<Int_t> staves_less8;
	vector<Int_t> staves_less4;
	vector<Int_t> staves_less3;
	vector<Int_t> staves_more3;
	vector<Int_t> staves_emulated_pp;
	vector<Int_t> staves_no_emulated_pp;
	vector<Int_t> staves_emulated_PbPb;
	vector<Int_t> staves_no_emulated_PbPb;
	
	//definition of vector<int> to store the different behaviour reduced
	vector<Int_t> reduced_staves_inactive;
	vector<Int_t> reduced_staves_less8;
	vector<Int_t> reduced_staves_less4;
	vector<Int_t> reduced_staves_less3;
	vector<Int_t> reduced_staves_more3;
	vector<Int_t> reduced_staves_emulated_pp;
	vector<Int_t> reduced_staves_no_emulated_pp;
	vector<Int_t> reduced_staves_emulated_PbPb;
	vector<Int_t> reduced_staves_no_emulated_PbPb;
		
	//definition of factor to determine threshold of majority of staves, i.e. fixing half of the staves
	Int_t threshold = 2;
	
	//definiton of other used variables 
	Float_t content;
	Int_t NA_counter = 0;
	//setup counter and vector of stave for each run
	for(Int_t j = 0; j < n_events; j++)
	{
		string check_cumulative = cumulative_label[j];
		Bool_t check = kFALSE;
		Int_t label_event_size = label_event.size();
		for(Int_t i = 0; i < label_event_size; i++)
		{
			string check_label = label_event[i];
			if (check_cumulative==check_label) check = kTRUE;
		}
		if(check)
		{
			for(Int_t k=0; k < n_stave; k++)
			{
				content = staves[k]->GetPointY(j-NA_counter);	
				if (content == 0) 
				{
					staves_inactive.push_back(k);
					counter_inactive++;
				}
				if (content >=1.0e-3) 
				{
					staves_more3.push_back(k);
					counter_more3++;
				}
				if (content >=1.0e-4 && content < 1.0e-3) 
				{
					staves_less3.push_back(k);
					counter_less3++;
				}
				if (content >=1.0e-8 && content < 1.0e-4) 
				{
					staves_less4.push_back(k);
					counter_less4++;
				}
				if (content < 1.0e-8 && content!=0) 
				{
					staves_less8.push_back(k);
					counter_less8++; 
				}
				if(layer == 0)
				{
				    if(content >= 1.9e-5 && content <= 2.1e-5)
				    {
				        staves_emulated_pp.push_back(k);
				        counter_emulated_pp++;
				    }
				    else staves_no_emulated_pp.push_back(k);
				    
				    if(content >= 6.6e-4 && content <= 7.2e-4)
				    {
				        staves_emulated_PbPb.push_back(k);
				        counter_emulated_PbPb++;
				    }
				    else staves_no_emulated_PbPb.push_back(k);
				}
				if(layer == 1)
				{
				    if(content >= 1.1e-5 && content <= 1.3e-5)
				    {
				        staves_emulated_pp.push_back(k);
				        counter_emulated_pp++;
				    }
				    else staves_no_emulated_pp.push_back(k);
				    
				    if(content >= 4.3e-4 && content <= 4.7e-4)
				    {
				        staves_emulated_PbPb.push_back(k);
				        counter_emulated_PbPb++;
				    }
				    else staves_no_emulated_PbPb.push_back(k);
				}
			    if(layer == 2)
				{
				    if(content >= 8.4e-6 && content <= 8.6e-6)
				    {
				        staves_emulated_pp.push_back(k);
				        counter_emulated_pp++;
				    }
				    else staves_no_emulated_pp.push_back(k);
				    
				    if(content >= 2.7e-4 && content <= 3.1e-4)
				    {
				        staves_emulated_PbPb.push_back(k);
				        counter_emulated_PbPb++;
				    }
				    else staves_no_emulated_PbPb.push_back(k);
				}
				if(layer == 3)
				{
				    if(content >= 5.2e-7 && content <= 5.4e-7)
				    {
				        staves_emulated_pp.push_back(k);
				        counter_emulated_pp++;
				    }
				    else staves_no_emulated_pp.push_back(k);
				    
				    if(content >= 8.7e-6 && content <= 9.6e-6)
				    {
				        staves_emulated_PbPb.push_back(k);
				        counter_emulated_PbPb++;
				    }
				    else staves_no_emulated_PbPb.push_back(k);
				}
				if(layer == 4)
				{
				    if(content >= 3.9e-7 && content <= 4.1e-7)
				    {
				        staves_emulated_pp.push_back(k);
				        counter_emulated_pp++;
				    }
				    else staves_no_emulated_pp.push_back(k);
				    
				    if(content >= 5.9e-6 && content <= 6.8e-6)
				    {
				        staves_emulated_PbPb.push_back(k);
				        counter_emulated_PbPb++;
				    }
				    else staves_no_emulated_PbPb.push_back(k);
				}
				if(layer == 5)
				{
				    if(content >= 2.5e-7 && content <= 2.7e-7)
				    {
				        staves_emulated_pp.push_back(k);
				        counter_emulated_pp++;
				    }
				    else staves_no_emulated_pp.push_back(k);
				    
				    if(content >= 3.3e-6 && content <= 3.8e-6)
				    {
				        staves_emulated_PbPb.push_back(k);
				        counter_emulated_PbPb++;
				    }
				    else staves_no_emulated_PbPb.push_back(k);
				}
				if(layer == 6)
				{
				    if(content >= 2.1e-7 && content <= 2.3e-7)
				    {
				        staves_emulated_pp.push_back(k);
				        counter_emulated_pp++;
				    }
				    else staves_no_emulated_pp.push_back(k);
				    
				    if(content >= 2.7e-6 && content <= 3.2e-6)
				    {
				        staves_emulated_PbPb.push_back(k);
				        counter_emulated_PbPb++;
				    }
				    else staves_no_emulated_PbPb.push_back(k);
				}
			}
								
			//cout on terminal about the analysis 
			//first check: if all the staves have FHR in the same range - cout a phrase
			if(counter_inactive==n_stave) { cout << "\e[31mrun " << cumulative_label[j] << "\t expert \t bad \t "<< layer_name << ", all stave: no response; \e[39m" << endl;}
			if(counter_less8==n_stave) {cout << "\e[32mrun " << cumulative_label[j] << "\t cosmic \t good \t " << layer_name << ", all staves: FHR < 10^-8; \e[39m" << endl;}		
			if(counter_less4==n_stave && counter_emulated_pp < n_stave/threshold && counter_emulated_PbPb < n_stave/threshold) {cout << "\e[34mrun " << cumulative_label[j] << "\t pp coll \t good \t " << layer_name << ", all staves: 10^-8 < FHR < 10^-4; \e[39m" <<endl;}
			if(counter_less3==n_stave && counter_emulated_pp < n_stave/threshold && counter_emulated_PbPb < n_stave/threshold) {cout << "\e[36mrun " << cumulative_label[j] << "\t Pb-Pb coll \t good \t " << layer_name << ", all staves: 10^-4 < FHR < 10^-3; \e[39m" << endl;}
			if(counter_more3==n_stave) {cout << "\e[30mrun " << cumulative_label[j] << "\t Pb-Pb coll \t bad \t " << layer_name << ", all staves: FHR > 10^-3; \e[39m" << endl;}
			if(counter_emulated_pp==n_stave) {cout << "\e[44mrun " << cumulative_label[j] << "\t pp emul \t good \t " << layer_name << ", all staves: in layer value; \e[0m" << endl;}
			if(counter_emulated_PbPb==n_stave) {cout << "\e[46mrun " << cumulative_label[j] << "\t Pb-Pb emul \t good \t " << layer_name << ", all staves: in layer value; \e[0m" << endl;}
			
			//second check: if the majority of the staves have a FHR value - cout a phrase with the staves with different behaviour
			//majority: inactive
			if(counter_inactive>=n_stave/threshold && counter_inactive!=n_stave) 
			{
				cout << "\e[31mrun " << cumulative_label[j] << " \t expert \t bad \t " << layer_name << ", stave ";
				if(staves_inactive.size()!=0)
				{
					extended_to_reduced_terminal(staves_inactive, reduced_staves_inactive);
					cout << ": no response";
				}
				if(staves_more3.size()!=0)
				{
					cout << ", stave ";
					extended_to_reduced_terminal(staves_more3, reduced_staves_more3);
					cout << ": FHR > 10^-3";
				}
				if(staves_less3.size()!=0)
				{
					cout << ", stave ";
					extended_to_reduced_terminal(staves_less3, reduced_staves_less3);
					cout << ": 10^-4 < FHR < 10^-3";
				}
				if(staves_less4.size()!=0)
				{
					cout << ", stave ";
					extended_to_reduced_terminal(staves_less4, reduced_staves_less4);
					cout << ": 10^-8 < FHR < 10^-4";
				}
				if(staves_less8.size()!=0)
				{
				    cout << ", stave ";
				    extended_to_reduced_terminal(staves_less8, reduced_staves_less8);
				    cout << ": FHR < 10^-8";
				}
				cout << "; \e[39m " << endl;
			}
			//majority: more3
			if(counter_more3>=n_stave/threshold && counter_more3!=n_stave && counter_inactive<n_stave/threshold) 
			{
				cout << "\e[30mrun " << cumulative_label[j] << " \t Pb-Pb col \t bad \t " << layer_name << ", stave ";
				if(staves_more3.size()!=0)
				{
					extended_to_reduced_terminal(staves_more3, reduced_staves_more3);
					cout << ": FHR > 10^-3";
				}
				if(staves_inactive.size()!=0)
				{
					cout << ", stave ";
					extended_to_reduced_terminal(staves_inactive, reduced_staves_inactive);
					cout << ": no response";
				}
				if(staves_less3.size()!=0)
				{
					cout << ", stave ";
					extended_to_reduced_terminal(staves_less3, reduced_staves_less3);
					cout << ": 10^-4 < FHR < 10^-3";
				}
				if(staves_less4.size()!=0)
				{
					cout << ", stave ";
					extended_to_reduced_terminal(staves_less4, reduced_staves_less4);
					cout << ": 10^-8 < FHR < 10^-4";
				}
				if(staves_less8.size()!=0)
				{
				    cout << ", stave ";
				    extended_to_reduced_terminal(staves_less8, reduced_staves_less8);
				    cout << ": FHR < 10^-8";
				}
				cout << "; \e[39m " << endl;
			}
					
			//majority: less3 - Pb-Pb candidate
			if(counter_less3>=n_stave/threshold && counter_less3!=n_stave && counter_inactive<n_stave/threshold && counter_more3<n_stave/threshold && counter_emulated_pp == 0 && counter_emulated_PbPb == 0) 
			{
				cout << "\e[36mrun " << cumulative_label[j] << "\t Pb-Pb col \t good \t " << layer_name << ", stave ";
				if(staves_less3.size()!=0)
				{
					extended_to_reduced_terminal(staves_less3, reduced_staves_less3);
					cout << ": 10^-4 < FHR < 10^-3";
				}
				if(staves_inactive.size()!=0)
				{
					cout << ", stave ";
					extended_to_reduced_terminal(staves_inactive, reduced_staves_inactive);
					cout << ": no response";
				}
				if(staves_more3.size()!=0)
				{
					cout << ", stave ";
					extended_to_reduced_terminal(staves_more3, reduced_staves_more3);
					cout << ": FHR > 10^-3";
				}
				if(staves_less4.size()!=0)
				{
					cout << ", stave ";
					extended_to_reduced_terminal(staves_less4, reduced_staves_less4);
					cout << ": 10^-8 < FHR < 10^-4";
				}
				if(staves_less8.size()!=0)
				{
				    cout << ", stave ";
				    extended_to_reduced_terminal(staves_less8, reduced_staves_less8);
				    cout << ": FHR < 10^-8";
				}
				cout << "; \e[39m " << endl;
			}
			
			//majority: less4, pp candidate
			if(counter_less4>=n_stave/threshold && counter_less4!=n_stave && counter_inactive<n_stave/threshold && counter_more3<n_stave/threshold && counter_less3<n_stave/threshold && counter_emulated_pp == 0 && counter_emulated_PbPb == 0) 
			{
				cout << "\e[34mrun " << cumulative_label[j] << "\t pp col \t good \t " << layer_name << ", stave ";
				if(staves_less4.size()!=0)
				{
					extended_to_reduced_terminal(staves_less4, reduced_staves_less4);
					cout << ": 10^-8 < FHR < 10^-4";
				}
				if(staves_inactive.size()!=0)
				{
					cout << ", stave ";
					extended_to_reduced_terminal(staves_inactive, reduced_staves_inactive);
					cout << ": no response";
				}
				if(staves_more3.size()!=0)
				{
					cout << ", stave ";
					extended_to_reduced_terminal(staves_more3, reduced_staves_more3);
					cout << ": FHR > 10^-3";
				}
				if(staves_less3.size()!=0)
				{
					cout << ", stave ";
					extended_to_reduced_terminal(staves_less3, reduced_staves_less3);
					cout << ": 10^-4 < FHR < 10^-3";
				}
				if(staves_less8.size()!=0)
				{
				    cout << ", stave ";
				    extended_to_reduced_terminal(staves_less8, reduced_staves_less8);
				    cout << ": FHR < 10^-8";
				}
				cout << "; \e[39m " << endl;
			}
			
			//majority: less8, cosmic good
			if(counter_less8>=n_stave/threshold && counter_less8!=n_stave && counter_inactive<n_stave/threshold &&  counter_more3<n_stave/threshold && counter_less3<n_stave/threshold && counter_less4<n_stave/threshold)  
			{
				cout << "\e[32mrun " << cumulative_label[j] << "\t cosmic \t good \t " << layer_name << ", stave ";
				if(staves_less8.size()!=0)
				{
				    extended_to_reduced_terminal(staves_less8, reduced_staves_less8);
				    cout << ": FHR < 10^-8";
				}
				if(staves_inactive.size()!=0)
				{
					cout << ", stave ";
					extended_to_reduced_terminal(staves_inactive, reduced_staves_inactive);
					cout << ": no response";
				}
				if(staves_more3.size()!=0)
				{
					cout << ", stave ";
					extended_to_reduced_terminal(staves_more3, reduced_staves_more3);
					cout << ": FHR > 10^-3";
				}
				if(staves_less3.size()!=0)
				{
					cout << ", stave ";
					extended_to_reduced_terminal(staves_less3, reduced_staves_less3);
					cout << ": 10^-4 < FHR < 10^-3";
				}
				if(staves_less4.size()!=0)
				{
					cout << ", stave ";
					extended_to_reduced_terminal(staves_less4, reduced_staves_less4);
					cout << ": 10^-8 < FHR < 10^-4";
				}
				cout << "; \e[39m " << endl;
			}
			//majority of the staves: emulated pp
			if(counter_emulated_pp >= n_stave/threshold && counter_emulated_pp != n_stave)
			{
			    cout << "\e[44mrun " << cumulative_label[j] << "\t pp emul \t good \t " << layer_name;
			    if(staves_emulated_pp.size()!=0)
			    {
			        cout << ", stave ";
			        extended_to_reduced_terminal(staves_emulated_pp, reduced_staves_emulated_pp);
			        cout << ": in layer values";
			    }
			    if(staves_no_emulated_pp.size()!=0)
			    {
			        cout << ", stave ";
			        extended_to_reduced_terminal(staves_no_emulated_pp, reduced_staves_no_emulated_pp);
			        cout << ": outside layer values";
			    }
			    cout << "; \e[0m " << endl;
			}
			//majority of staves: emulated PbPb
			if(counter_emulated_PbPb >= n_stave/threshold && counter_emulated_PbPb != n_stave)
			{
			    cout << "\e[46mrun " << cumulative_label[j] << "\t Pb-Pb emul \t good \t " << layer_name;
			    if(staves_emulated_PbPb.size()!=0)
			    {
			        cout << ", stave ";
			        extended_to_reduced_terminal(staves_emulated_PbPb, reduced_staves_emulated_PbPb);
			        cout << ": in layer values";
			    }
			    if(staves_no_emulated_PbPb.size()!=0)
			    {
			        cout << ", stave ";
			        extended_to_reduced_terminal(staves_no_emulated_PbPb, reduced_staves_no_emulated_PbPb);
			        cout << ": outside layer values";
			    }
			    cout << "; \e[0m " << endl;
			}
			//other cases
			if(counter_inactive < n_stave/threshold && counter_more3 < n_stave/threshold && counter_less3 < n_stave/threshold && counter_less4 < n_stave/threshold && counter_less8 < n_stave/threshold)
			{
				cout << "\e[30mrun " << cumulative_label[j] << " \t expert \t exp \t " << layer_name;
				if(staves_inactive.size()!=0)
				{
					cout << ", stave ";
					extended_to_reduced_terminal(staves_inactive, reduced_staves_inactive);
					cout << ": no response";
				}
				if(staves_more3.size()!=0)
				{
					cout << ", stave ";
					extended_to_reduced_terminal(staves_more3, reduced_staves_more3);
					cout << ": FHR > 10^-3";
				}
				if(staves_less3.size()!=0)
				{
					cout << ", stave ";
					extended_to_reduced_terminal(staves_less3, reduced_staves_less3);
					cout << ": 10^-4 < FHR < 10^-3";
				}
				if(staves_less4.size()!=0)
				{
					cout << ", stave ";
					extended_to_reduced_terminal(staves_less4, reduced_staves_less4);
					cout << ": 10^-8 < FHR < 10^-4";
				}
				if(staves_less8.size()!=0)
				{
					cout << ", stave ";
					extended_to_reduced_terminal(staves_less8, reduced_staves_less8);
					cout << ": FHR < 10^-8";
				}
				cout << "; \e[39m " << endl;
			}
		}
		
		else 
		{
			cout <<"\e[35mrun " << cumulative_label[j] << "\t expert \t N/A \t \t \e[39m" << endl;
			NA_counter++;
		}	
		
		//put all the counter to zero for a new run
		counter_less8=0;
		counter_less4=0;
		counter_less3=0;
		counter_more3=0;
		counter_inactive=0;
		counter_emulated_pp=0;
		counter_emulated_PbPb=0;
		
		//reset vector of staves for a new run
		staves_inactive.clear();
		staves_more3.clear();
		staves_less3.clear();
		staves_less4.clear();
		staves_less8.clear();
		staves_emulated_pp.clear();
		staves_no_emulated_pp.clear();
		staves_emulated_PbPb.clear();
		staves_no_emulated_PbPb.clear();
		
		reduced_staves_inactive.clear();
		reduced_staves_more3.clear();
		reduced_staves_less3.clear();
		reduced_staves_less4.clear();
		reduced_staves_less8.clear();
		reduced_staves_emulated_pp.clear();
		reduced_staves_no_emulated_pp.clear();
		reduced_staves_emulated_PbPb.clear();
		reduced_staves_no_emulated_PbPb.clear();
	} 
};

//analysis for the repository
void create_repository(Int_t layer, Int_t n_stave, Int_t n_events, vector<string> cumulative_label, vector<string> label_event, vector<TGraph*> staves, const char* layer_name, string final_comment, vector<string>& comments, vector<Int_t>& build_layer_status, vector<Int_t>& build_layer_collision)
{

	
	//definition of counters to help the division
	Int_t counter_inactive = 0;
	Int_t counter_less8 = 0; //counter for FHR < 10^-8
	Int_t counter_less4 = 0; //counter for 10^-8 < FHR < 10^-4 
	Int_t counter_less3 = 0; //counter for 10^-4 < FHR < 10^-3
	Int_t counter_more3 = 0; //counter for FHR > 10^-3
	Int_t counter_emulated_pp = 0; //counter for emulated pp
	Int_t counter_emulated_PbPb = 0; //counter for emulated PbPb
	
	//definition of vector<int> to store the different behaviour extended
	vector<Int_t> staves_inactive;
	vector<Int_t> staves_less8;
	vector<Int_t> staves_less4;
	vector<Int_t> staves_less3;
	vector<Int_t> staves_more3;
	vector<Int_t> staves_emulated_pp;
	vector<Int_t> staves_no_emulated_pp;
	vector<Int_t> staves_emulated_PbPb;
	vector<Int_t> staves_no_emulated_PbPb;
	
	//definition of vector<int> to store the different behaviour reduced
	vector<Int_t> reduced_staves_inactive;
	vector<Int_t> reduced_staves_less8;
	vector<Int_t> reduced_staves_less4;
	vector<Int_t> reduced_staves_less3;
	vector<Int_t> reduced_staves_more3;
    vector<Int_t> reduced_staves_emulated_pp;
    vector<Int_t> reduced_staves_no_emulated_pp;
    vector<Int_t> reduced_staves_emulated_PbPb;   
    vector<Int_t> reduced_staves_no_emulated_PbPb;
    		
	//definition of factor to determine threshold of majority of staves, i.e. fixing half of the staves
	Int_t threshold = 2;
	
	//definiton of other used variables 
	Float_t content;
	Int_t NA_counter = 0;
	//setup counter and vector of stavef for each run
	for(Int_t j = 0; j < n_events; j++)
	{
		string check_cumulative = cumulative_label[j];
		Bool_t check = kFALSE;
		Int_t label_event_size = label_event.size();
		for(Int_t i = 0; i < label_event_size; i++)
		{
			string check_label = label_event[i];
			if (check_cumulative==check_label) check = kTRUE;
		}
		if(check)
		{
			for(Int_t k=0; k < n_stave; k++)
			{
				content = staves[k]->GetPointY(j-NA_counter);	
				if (content == 0) 
				{
					staves_inactive.push_back(k);
					counter_inactive++;
				}
				if (content >=1.0e-3) 
				{
					staves_more3.push_back(k);
					counter_more3++;
				}
				if (content >=1.0e-4 && content < 1.0e-3) 
				{
					staves_less3.push_back(k);
					counter_less3++;
				}
				if (content >=1.0e-8 && content < 1.0e-4) 
				{
					staves_less4.push_back(k);
					counter_less4++;
				}
				if (content < 1.0e-8 && content!=0) 
				{
					staves_less8.push_back(k);
					counter_less8++; 
				}
				if(layer == 0)
				{
				    if(content >= 1.9e-5 && content <= 2.1e-5)
				    {
				        staves_emulated_pp.push_back(k);
				        counter_emulated_pp++;
				    }
				    else staves_no_emulated_pp.push_back(k);
				    
				    if(content >= 6.6e-4 && content <= 7.2e-4)
				    {
				        staves_emulated_PbPb.push_back(k);
				        counter_emulated_PbPb++;
				    }
				    else staves_no_emulated_PbPb.push_back(k);
				}
				if(layer == 1)
				{
				    if(content >= 1.1e-5 && content <= 1.3e-5)
				    {
				        staves_emulated_pp.push_back(k);
				        counter_emulated_pp++;
				    }
				    else staves_no_emulated_pp.push_back(k);
				    
				    if(content >= 4.3e-4 && content <= 4.7e-4)
				    {
				        staves_emulated_PbPb.push_back(k);
				        counter_emulated_PbPb++;
				    }
				    else staves_no_emulated_PbPb.push_back(k);
				}
			    if(layer == 2)
				{
				    if(content >= 8.4e-6 && content <= 8.6e-6)
				    {
				        staves_emulated_pp.push_back(k);
				        counter_emulated_pp++;
				    }
				    else staves_no_emulated_pp.push_back(k);
				    
				    if(content >= 2.7e-4 && content <= 3.1e-4)
				    {
				        staves_emulated_PbPb.push_back(k);
				        counter_emulated_PbPb++;
				    }
				    else staves_no_emulated_PbPb.push_back(k);
				}
				if(layer == 3)
				{
				    if(content >= 5.2e-7 && content <= 5.4e-7)
				    {
				        staves_emulated_pp.push_back(k);
				        counter_emulated_pp++;
				    }
				    else staves_no_emulated_pp.push_back(k);
				    
				    if(content >= 8.7e-6 && content <= 9.6e-6)
				    {
				        staves_emulated_PbPb.push_back(k);
				        counter_emulated_PbPb++;
				    }
				    else staves_no_emulated_PbPb.push_back(k);
				}
				if(layer == 4)
				{
				    if(content >= 3.9e-7 && content <= 4.1e-7)
				    {
				        staves_emulated_pp.push_back(k);
				        counter_emulated_pp++;
				    }
				    else staves_no_emulated_pp.push_back(k);
				    
				    if(content >= 5.9e-6 && content <= 6.8e-6)
				    {
				        staves_emulated_PbPb.push_back(k);
				        counter_emulated_PbPb++;
				    }
				    else staves_no_emulated_PbPb.push_back(k);
				}
				if(layer == 5)
				{
				    if(content >= 2.5e-7 && content <= 2.7e-7)
				    {
				        staves_emulated_pp.push_back(k);
				        counter_emulated_pp++;
				    }
				    else staves_no_emulated_pp.push_back(k);
				    
				    if(content >= 3.3e-6 && content <= 3.8e-6)
				    {
				        staves_emulated_PbPb.push_back(k);
				        counter_emulated_PbPb++;
				    }
				    else staves_no_emulated_PbPb.push_back(k);
				}
				if(layer == 6)
				{
				    if(content >= 2.1e-7 && content <= 2.3e-7)
				    {
				        staves_emulated_pp.push_back(k);
				        counter_emulated_pp++;
				    }
				    else staves_no_emulated_pp.push_back(k);
				    
				    if(content >= 2.7e-6 && content <= 3.2e-6)
				    {
				        staves_emulated_PbPb.push_back(k);
				        counter_emulated_PbPb++;
				    }
				    else staves_no_emulated_PbPb.push_back(k);
				}
			}
			   
			//layer status: 0 = expert check 1 = bad; 2 = good; 3 = N/A
			//layer collision: 0 = expert check; 1 = cosmic; 2 = pp coll; 3 = Pb-Pb coll; 4 = pp emul; 5: Pb-Pb emul;
			//create comments for the repository
			//first check: if all the staves have FHR in the same range
			if(counter_inactive==n_stave) 
			{
				final_comment.append(layer_name);
				final_comment.append(", all stave: no response; ");
				build_layer_status.push_back(1);
				build_layer_collision.push_back(0);
			}
			if(counter_less8==n_stave) 
			{
				final_comment.append(layer_name);
				final_comment.append(", all staves: FHR < 10^-8; ");
				build_layer_status.push_back(2);
				build_layer_collision.push_back(1);	
			}		
			
			if(counter_less4==n_stave && counter_emulated_pp < n_stave/threshold && counter_emulated_PbPb < n_stave/threshold) 
			{
				final_comment.append(layer_name);
				final_comment.append(", all staves: 10^-8 < FHR < 10^-4; ");
				build_layer_status.push_back(2);
				build_layer_collision.push_back(2);
			}
			if(counter_less3==n_stave && counter_emulated_pp < n_stave/threshold && counter_emulated_PbPb < n_stave/threshold) 
			{
				final_comment.append(layer_name);
				final_comment.append(", all staves: 10^-4 < FHR < 10^-3; ");
				build_layer_status.push_back(2);
				build_layer_collision.push_back(3);
			}
			if(counter_more3==n_stave) 
			{
				final_comment.append(layer_name);
				final_comment.append(", all staves: FHR > 10^-3; ");
				build_layer_status.push_back(1);
				build_layer_collision.push_back(0);
			}
			if(counter_emulated_pp==n_stave)
			{
			    final_comment.append(layer_name);
			    final_comment.append(", all staves: inside layer values; ");
			    build_layer_status.push_back(2);
			    build_layer_collision.push_back(4);
			}
			if(counter_emulated_PbPb==n_stave)
			{
			    final_comment.append(layer_name);
			    final_comment.append(", all staves: inside layer values; ");
			    build_layer_status.push_back(2);
			    build_layer_collision.push_back(5);
			}	
			//second check: if the majority of the staves have a FHR value - create a phrase with the staves with different behaviour
			//majority: inactive
			if(counter_inactive>=n_stave/threshold && counter_inactive!=n_stave) 
			{
				final_comment.append(layer_name);
				final_comment.append(", stave ");
				build_layer_status.push_back(1);
				build_layer_collision.push_back(0);
				if(staves_inactive.size()!=0)
				{
					extended_to_reduced_repository(staves_inactive, reduced_staves_inactive, final_comment);
					final_comment.append(": no response");
				}
				if(staves_more3.size()!=0)
				{
					final_comment.append(", stave ");
					extended_to_reduced_repository(staves_more3, reduced_staves_more3, final_comment);
					final_comment.append(": FHR > 10^-3");
				}
				if(staves_less3.size()!=0)
				{
					final_comment.append(", stave ");
					extended_to_reduced_repository(staves_less3, reduced_staves_less3, final_comment);
					final_comment.append(": 10^-4 < FHR < 10^-3");
				}
				if(staves_less4.size()!=0)
				{
					final_comment.append(", stave ");
					extended_to_reduced_repository(staves_less4, reduced_staves_less4, final_comment);
					final_comment.append(": 10^-8 < FHR < 10^-4");
				}
				if(staves_less8.size()!=0)
				{
				    final_comment.append(", stave ");
				    extended_to_reduced_repository(staves_less8, reduced_staves_less8, final_comment);
				    final_comment.append(": FHR < 10^-8");
				}				
				final_comment.append("; ");
			}
			//majority: more3
			if(counter_more3>=n_stave/threshold && counter_more3!=n_stave && counter_inactive<n_stave/threshold) 
			{
				final_comment.append(layer_name);
				final_comment.append(", stave ");
				build_layer_status.push_back(1);
				build_layer_collision.push_back(0);
				if(staves_more3.size()!=0)
				{
					extended_to_reduced_repository(staves_more3, reduced_staves_more3, final_comment);
					final_comment.append(": FHR > 10^-3");
				}
				if(staves_inactive.size()!=0)
				{
					final_comment.append(", stave ");
					extended_to_reduced_repository(staves_inactive, reduced_staves_inactive, final_comment);
					final_comment.append(": no response");
				}
				if(staves_less3.size()!=0)
				{
					final_comment.append(", stave ");
					extended_to_reduced_repository(staves_less3, reduced_staves_less3, final_comment);
					final_comment.append(": 10^-4 < FHR < 10^-3");
				}
				if(staves_less4.size()!=0)
				{
					final_comment.append(", stave ");
					extended_to_reduced_repository(staves_less4, reduced_staves_less4, final_comment);
					final_comment.append(": 10^-8 < FHR < 10^-4");
				}
				if(staves_less8.size()!=0)
				{
				    final_comment.append(", stave ");
				    extended_to_reduced_repository(staves_less8, reduced_staves_less8, final_comment);
				    final_comment.append(": FHR < 10^-8");
				}
				final_comment.append("; ");
			}
			//majority: less3 - Pb-Pb candidate
			if(counter_less3>=n_stave/threshold && counter_less3!=n_stave && counter_inactive<n_stave/threshold && counter_more3<n_stave/threshold && counter_emulated_pp < n_stave/threshold && counter_emulated_PbPb < n_stave/threshold) 
			{
				final_comment.append(layer_name);
				final_comment.append(", stave ");
				build_layer_status.push_back(2);
				build_layer_collision.push_back(3);
				if(staves_less3.size()!=0)
				{
					extended_to_reduced_repository(staves_less3, reduced_staves_less3, final_comment);
					final_comment.append(": 10^-4 < FHR < 10^-3");
				}
				if(staves_inactive.size()!=0)
				{
					final_comment.append(", stave ");
					extended_to_reduced_repository(staves_inactive, reduced_staves_inactive, final_comment);
					final_comment.append(": no response");
				}
				if(staves_more3.size()!=0)
				{
					final_comment.append(", stave ");
					extended_to_reduced_repository(staves_more3, reduced_staves_more3, final_comment);
					final_comment.append(": FHR > 10^-3");
				}
				if(staves_less4.size()!=0)
				{
					final_comment.append(", stave ");
					extended_to_reduced_repository(staves_less4, reduced_staves_less4, final_comment);
					final_comment.append(": 10^-8 < FHR < 10^-4");
				}
				if(staves_less8.size()!=0)
				{
				    final_comment.append(", stave ");
				    extended_to_reduced_repository(staves_less8, reduced_staves_less8, final_comment);
				    final_comment.append(": FHR < 10^-8");
				}
			}
			//majority: less4, pp candidate
			if(counter_less4>=n_stave/threshold && counter_less4!=n_stave && counter_inactive<n_stave/threshold && counter_more3<n_stave/threshold && counter_less3<n_stave/threshold && counter_emulated_pp < n_stave/threshold && counter_emulated_PbPb < n_stave/threshold) 
			{
				final_comment.append(layer_name);
				final_comment.append(", stave ");
				build_layer_status.push_back(2);
				build_layer_collision.push_back(2);
				if(staves_less4.size()!=0)
				{
					extended_to_reduced_repository(staves_less4, reduced_staves_less4, final_comment);
					final_comment.append(": 10^-8 < FHR < 10^-4");
				}
				if(staves_inactive.size()!=0)
				{
					final_comment.append(", stave ");
					extended_to_reduced_repository(staves_inactive, reduced_staves_inactive, final_comment);
					final_comment.append(": no response");
				}
				if(staves_more3.size()!=0)
				{
					final_comment.append(", stave ");
					extended_to_reduced_repository(staves_more3, reduced_staves_more3, final_comment);
					final_comment.append(": FHR > 10^-3");
				}
				if(staves_less3.size()!=0)
				{
					final_comment.append(", stave ");
					extended_to_reduced_repository(staves_less3, reduced_staves_less3, final_comment);
					final_comment.append(": 10^-4 < FHR < 10^-3");
				}
				if(staves_less8.size()!=0)
				{
				    final_comment.append(", stave ");
				    extended_to_reduced_repository(staves_less8, reduced_staves_less8, final_comment);
				    final_comment.append(": FHR < 10^-8");
				}
				final_comment.append("; ");
			}
			
			//majority: less8, cosmic good
			if(counter_less8>=n_stave/threshold && counter_less8!=n_stave && counter_inactive<n_stave/threshold && counter_more3<n_stave/threshold && counter_less3<n_stave/threshold && counter_less4<n_stave/threshold)  
			{
				final_comment.append(layer_name);
				final_comment.append(", stave ");
				build_layer_status.push_back(2);
				build_layer_collision.push_back(1);
				if(staves_less8.size()!=0)
				{
				    extended_to_reduced_repository(staves_less8, reduced_staves_less8, final_comment);
				    final_comment.append(": FHR < 10^-8");
				}
				if(staves_inactive.size()!=0)
				{
					final_comment.append(", stave ");
					extended_to_reduced_repository(staves_inactive, reduced_staves_inactive, final_comment);
					final_comment.append(": no response");
				}
				if(staves_more3.size()!=0)
				{
					final_comment.append(", stave ");
					extended_to_reduced_repository(staves_more3, reduced_staves_more3, final_comment);
					final_comment.append(": FHR > 10^-3");
				}
				if(staves_less3.size()!=0)
				{
					final_comment.append(", stave "); 
					extended_to_reduced_repository(staves_less3, reduced_staves_less3, final_comment);
					final_comment.append(": 10^-4 < FHR < 10^-3");
				}
				if(staves_less4.size()!=0)
				{
					final_comment.append(", stave ");
					extended_to_reduced_repository(staves_less4, reduced_staves_less4, final_comment);
					final_comment.append(": 10^-8 < FHR < 10^-4");
				}
				final_comment.append("; ");
			}
			//majority: pp emulated
			if(counter_emulated_pp >= n_stave/threshold && counter_emulated_pp!= n_stave)
			{
			    final_comment.append(layer_name);
				final_comment.append(", stave ");
				build_layer_status.push_back(2);
				build_layer_collision.push_back(4);
				if(staves_emulated_pp.size()!=0)
				{
				    extended_to_reduced_repository(staves_emulated_pp, reduced_staves_emulated_pp, final_comment);
				    final_comment.append(": inside layer values");
				}
				if(staves_no_emulated_pp.size()!=0)
				{
				    final_comment.append(", stave ");
				    extended_to_reduced_repository(staves_no_emulated_pp, reduced_staves_no_emulated_pp, final_comment);
				    final_comment.append(": outside layer values");
				}
				final_comment.append("; ");
			}
			//majority: PbPb emulated
			if(counter_emulated_PbPb >= n_stave/threshold && counter_emulated_PbPb!= n_stave)
			{
			    final_comment.append(layer_name);
				final_comment.append(", stave ");
				build_layer_status.push_back(2);
				build_layer_collision.push_back(5);
				if(staves_emulated_PbPb.size()!=0)
				{
				    extended_to_reduced_repository(staves_emulated_PbPb, reduced_staves_emulated_PbPb, final_comment);
				    final_comment.append(": inside layer values");
				}
				if(staves_no_emulated_PbPb.size()!=0)
				{
				    final_comment.append(", stave ");
				    extended_to_reduced_repository(staves_no_emulated_PbPb, reduced_staves_no_emulated_PbPb, final_comment);
				    final_comment.append(": outside layer values");
				}
				final_comment.append("; ");
			}
			//other cases
			if(counter_inactive < n_stave/threshold && counter_more3 < n_stave/threshold && counter_less3 < n_stave/threshold && counter_less4 < n_stave/threshold && counter_less8 < n_stave/threshold)
			{
				final_comment.append(layer_name);
				build_layer_status.push_back(0);
				build_layer_collision.push_back(0);
				if(staves_inactive.size()!=0)
				{
					final_comment.append(", stave ");
					extended_to_reduced_repository(staves_inactive, reduced_staves_inactive, final_comment);
					final_comment.append(": no response");
				}
				if(staves_more3.size()!=0)
				{
					final_comment.append(", stave ");
					extended_to_reduced_repository(staves_more3, reduced_staves_more3, final_comment);
					final_comment.append(": FHR > 10^-3");
				}
				if(staves_less3.size()!=0)
				{
					final_comment.append(", stave "); 
					extended_to_reduced_repository(staves_less3, reduced_staves_less3, final_comment);
					final_comment.append(": 10^-4 < FHR < 10^-3");
				}
				if(staves_less4.size()!=0)
				{
					final_comment.append(", stave ");
					extended_to_reduced_repository(staves_less4, reduced_staves_less4, final_comment);
					final_comment.append(": 10^-8 < FHR < 10^-4");
				}
				if(staves_less8.size()!=0)
				{
					final_comment.append(", stave ");
					extended_to_reduced_repository(staves_less8, reduced_staves_less8, final_comment);
					final_comment.append(": FHR < 10^-8");
				}
				final_comment.append("; ");
			}
			
			comments.push_back(final_comment);
			final_comment = "";
		}
		
		else 
		{
			final_comment = "";
			comments.push_back(final_comment);
			build_layer_status.push_back(3);
			build_layer_collision.push_back(0);
			NA_counter++;
		}
		//put all the counter to zero for a new run
		counter_less8=0;
		counter_less4=0;
		counter_less3=0;
		counter_more3=0;
		counter_inactive=0;
		counter_emulated_pp=0;
		counter_emulated_PbPb=0;
		
		//reset vector of staves for a new run
		staves_inactive.clear();
		staves_more3.clear();
		staves_less3.clear();
		staves_less4.clear();
		staves_less8.clear();
		staves_emulated_pp.clear();
		staves_no_emulated_pp.clear();
		staves_emulated_PbPb.clear();
		staves_no_emulated_PbPb.clear();
		
		reduced_staves_inactive.clear();
		reduced_staves_more3.clear();
		reduced_staves_less3.clear();
		reduced_staves_less4.clear();
		reduced_staves_less8.clear();
		reduced_staves_emulated_pp.clear();
		reduced_staves_no_emulated_pp.clear();
		reduced_staves_emulated_PbPb.clear();
		reduced_staves_no_emulated_PbPb.clear();
	
	} 
};

//main scope
void CreateRepositoryFHR()
{
	gROOT->SetBatch(kTRUE);
	
	//definition of elements that will be used inside the loop
	TCanvas* c;
	TH1F* hist;
	TFile* input_file;
	vector<vector<string>> label_event;
	vector<string> label_prv;
	vector<string> cumulative_label;
	string label;
	Int_t n_events[7];
	Int_t number_staves[7] = {12, 16, 20, 24, 30, 42, 48};
	Int_t n_stave;
	const char* layer_name;
	
	//definition of variables to choose the file/section/layer of interest
	Int_t run_lower, run_upper; 
	Int_t sector_analysis;
	Int_t layer_lower, layer_upper;
	string layer_choice="";
	Int_t layer_number; 
	Bool_t terminal = kFALSE;
	string terminal_check = "";
	Bool_t repository = kFALSE;
	string repository_check = "";
	Bool_t check_file=kTRUE;
	
	
	//choose the file and sector to analyze, and choose if create the repository
	cout << "Insert the lower edge of run interval " << endl; 
	cin >> run_lower; 
	cout << "Insert the upper edge of run interval " << endl; 
	cin >> run_upper; 
	cout << "Do you want to check the results on terminal? (y/n) " << endl;
	cin >> terminal_check;
	if(terminal_check=="y"||terminal_check=="Y") 
	{
		terminal = kTRUE;
		cout << "The results will be print on terminal " << endl;
	}
	else if (terminal_check=="n"||terminal_check=="N")
	{
		cout << "The results will NOT be print on terminal " << endl;
	}
	else 
	{
		cout << "Wrong choise. Restart the analysis" << endl;
		exit(1);
	}
	cout << "Do you want to create the repository file for ALL the layers? (y/n) " << endl; 
	cin >> repository_check;
	if(repository_check =="Y" || repository_check=="y") 
	{
		cout << "The repository will be created " << endl;
		repository = kTRUE;
	}
	else if(repository_check=="n" || repository_check=="N") 
	{
		cout << "The repository will NOT be created" << endl;
	}
	else 
	{
		cout << "Wrong choise. Restart the analysis. " << endl;
		exit(1);
	}
	
	//build the cumulative run list
	for(Int_t i=0; i<7; i++)
	{
		if(i<3)
		{
			if(gSystem->AccessPathName(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_w_error_data.root", i, run_lower, run_upper))&&gSystem->AccessPathName(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d.root", i, run_lower, run_upper)))
			{
				cout << "The file for layer " << i << " does not exist. Check the run extremes " << endl; 
				exit(1);
			}
			else 
			{
				input_file = TFile::Open(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_w_error_data.root", i, run_lower, run_upper));
				if(!input_file) input_file = TFile::Open(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d.root", i, run_lower, run_upper));
				c = (TCanvas*)input_file->Get("canvas");
				if(!c) c = (TCanvas*)input_file->Get(Form("Layer%d_fakehitrate_w_error_and_trig_data", i));
				
				hist = (TH1F*) c->FindObject(Form("hfake_L%d", i));
				n_events[i] = hist->GetNbinsX();
				for(Int_t j=0; j < n_events[i]; j++)
				{
				 	label = hist->GetXaxis()->GetBinLabel(j+1);
				   	label.erase(0, 3);
				   	label_prv.push_back(label);
				   	cumulative_label.push_back(label);
				}
				label_event.push_back(label_prv);
				label_prv.clear();	
			}	
		}
		else
		{
			if(gSystem->AccessPathName(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_w_error_data_HSLower.root", i, run_lower, run_upper))&&gSystem->AccessPathName(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_HSLower.root", i, run_lower, run_upper)))
			{
                    cout << "The file for layer " << i << " HS lower does not exist. Check the run extremes. " << endl;
				    exit(1);
			}
			else
			{
				input_file = TFile::Open(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_w_error_data_HSLower.root", i, run_lower, run_upper));
				if(!input_file) input_file = TFile::Open(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_HSLower.root", i, run_lower, run_upper)); 
				c = (TCanvas*)input_file->Get("canvas_HSLower");
				if(!c) c = (TCanvas*)input_file->Get(Form("Layer%d_fakehitrate_w_error_and_trig_data_HSLower", i));
				
				hist = (TH1F*) c->FindObject(Form("hfake_L%d", i));
			    n_events[i] = hist->GetNbinsX();
			    for(Int_t j=0; j < n_events[i]; j++)
			    {
			    	label = hist->GetXaxis()->GetBinLabel(j+1);
			    	label.erase(0, 3);
			    	label_prv.push_back(label);
			    	cumulative_label.push_back(label);
			    }
			    label_event.push_back(label_prv);
			    label_prv.clear();
			}
		}	
	}
	
	sort(cumulative_label.begin(), cumulative_label.end());	
	cumulative_label.erase(unique(cumulative_label.begin(), cumulative_label.end()), cumulative_label.end());
	
	//output on terminal
	Int_t total_events = cumulative_label.size();
	if(terminal)
	{
		cout << "For terminal output: press 0 for IB analysis, 1 for OB analysis or 2 for both " << endl;
		cin >> sector_analysis;
		if(sector_analysis!=0 && sector_analysis!=1 && sector_analysis!=2)
		{
			cout << "Wrong choise. Restart the analysis" << endl;
			exit(1);
		} 
		if(sector_analysis==0) {layer_lower = 0; layer_upper = 3;} 
		if(sector_analysis==1) {layer_lower = 3; layer_upper = 7;}
		if(sector_analysis==2) {layer_lower = 0; layer_upper = 7;}
	
		if(sector_analysis==0 || sector_analysis==1)
		{
			cout << "Would you like to analyze a single layer? (y/n) " << endl;
			cin >> layer_choice; 
			if(layer_choice=="y" || layer_choice=="Y")
			{
				cout << "Insert the number of the layer you want to analyze " << endl;
				cin >> layer_number;
				if(sector_analysis==0 && layer_number!=0 && layer_number!=1 && layer_number!=2)
				{
					cout << "Wrong choise. Restart the analysis " << endl;
					exit(1);
				}
				else {layer_lower = layer_number; layer_upper = layer_number+1;}
			
				if(sector_analysis==1 && layer_number!=3 && layer_number!=4 && layer_number!=5 && layer_number!=6)
				{
					cout << "Wrong choise. Restart the analysis " << endl;
					exit(1);
				}
				else {layer_lower = layer_number; layer_upper = layer_number+1;}	
			}
			else if(layer_choice=="n" || layer_choice=="N") cout << "Proceding with the analysis af all layers in the sector selected " << endl;
			else {cout << "wrong choise. Restart the analysis " << endl; exit(1);}
		}
		for(Int_t i=layer_lower; i < layer_upper; i++)
		{
			if(i < 3)
			{
				if(gSystem->AccessPathName(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_w_error_data.root", i, run_lower, run_upper))&&gSystem->AccessPathName(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d.root", i, run_lower, run_upper)))
				{
				    cout << "The file for layer " << i << " does not exist. Check the run extremes " << endl; 
					exit(1);	
				}
				else 
				{
					input_file = TFile::Open(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_w_error_data.root", i, run_lower, run_upper));
					if(!input_file) input_file = TFile::Open(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d.root", i, run_lower, run_upper));
					c = (TCanvas*)input_file->Get("canvas");
					if(!c) c = (TCanvas*)input_file->Get(Form("Layer%d_fakehitrate_w_error_and_trig_data", i));
					
				    n_stave = number_staves[i];
				    vector<TGraph*> staves;
				    for(Int_t j=0; j < n_stave; j++) 
				    {
				    	TGraph* aux_stave = (TGraph*) c->GetPrimitive((Form("gr_L%d_stave%d_HS0", i, j)));
				    	staves.push_back(aux_stave);
				    }
				    layer_name = Form("L%d", i);
 				    terminal_output(i, n_stave, total_events, cumulative_label, label_event[i], staves, layer_name);
    				cout << "operation on layer " << i << " concluded " << endl;
    				staves.clear();
				}
			}
			else
			{
				if(gSystem->AccessPathName(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_w_error_data_HSLower.root", i, run_lower, run_upper))&&gSystem->AccessPathName(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_HSLower.root", i, run_lower, run_upper)))
				{
				    cout << "The file for layer " << i << " HS lower does not exist. Check the run extremes. " << endl;
					exit(1);
				}
				else
				{
					input_file = TFile::Open(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_w_error_data_HSLower.root", i, run_lower, run_upper));
					if(!input_file) input_file = TFile::Open(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_HSLower.root", i, run_lower, run_upper));
					c = (TCanvas*)input_file->Get("canvas_HSLower");
					if(!c) c = (TCanvas*)input_file->Get(Form("Layer%d_fakehitrate_w_error_and_trig_data_HSLower", i));
					
					n_stave = number_staves[i];
					vector<TGraph*> staves_lower;
				    for(Int_t j=0; j < n_stave; j++) 
				    {
			    		TGraph* aux_stave = (TGraph*) c->GetPrimitive((Form("gr_L%d_stave%d_HS0", i, j)));
			    		staves_lower.push_back(aux_stave);
		    		}
		    		layer_name = Form("L%dL", i);
		    		terminal_output(i, n_stave, total_events, cumulative_label, label_event[i], staves_lower, layer_name);
					cout << "operation on layer " << i << " lower concluded " << endl;
					staves_lower.clear();
					
				}
				
				if(gSystem->AccessPathName(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_w_error_data_HSUpper.root", i, run_lower, run_upper))&&gSystem->AccessPathName(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_HSUpper.root", i, run_lower, run_upper)))
				{
					cout << "The file for layer " << i << " HS upper does not exist. Check the run extremes. " << endl;
					exit(1);    
				}
				else
				{
					input_file = TFile::Open(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_w_error_data_HSUpper.root", i, run_lower, run_upper));
					if(!input_file) input_file = TFile::Open(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_HSUpper.root", i, run_lower, run_upper)); 
					c = (TCanvas*)input_file->Get("canvas_HSUpper");
					if(!c) c = (TCanvas*)input_file->Get(Form("Layer%d_fakehitrate_w_error_and_trig_data_HSUpper", i));
					
					n_stave = number_staves[i];
					vector<TGraph*> staves_upper;
					for(Int_t j=0; j < n_stave; j++) 
					{
					  TGraph* aux_stave = (TGraph*) c->GetPrimitive((Form("gr_L%d_stave%d_HS1", i, j)));
					  staves_upper.push_back(aux_stave);
				    }
				    layer_name = Form("L%dU", i);
				    terminal_output(i, n_stave, total_events, cumulative_label, label_event[i], staves_upper, layer_name);
				    cout << "operation on layer " << i << " upper concluded " << endl;
				    staves_upper.clear();
				}
			}
		}
	}
	
	//create repository
	if(repository)
	{
		//definition of variables used for the repository creation
		vector<string> build_vcomments;
		vector<vector<string>> comments;
		vector<Int_t> build_layer_status;
		vector<Int_t> build_layer_collision;
		vector<Int_t> global_status; 
		vector<Int_t> global_collision;
		vector<string> global_status_descriptive;
		vector<vector<string>> layer_status_descriptive;
		vector<vector<Int_t>> layer_status;
		vector<vector<Int_t>> layer_collision;
		vector<string> definitive_comments;
		string build_definitive_comments; 
		Int_t min_status;
		vector<string> build_layer_status_descriptive;
		Int_t collision_type;
		vector<string> global_collision_descriptive;
		
		//start analysis
		for(Int_t i=0; i < 7; i++)
		{
			string final_comment = "";
			if(i < 3)
			{
				if(gSystem->AccessPathName(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_w_error_data.root", i, run_lower, run_upper))&&gSystem->AccessPathName(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d.root", i, run_lower, run_upper)))
				{
				    cout << "The file for layer " << i << " does not exist. Check the run extremes " << endl; 
					exit(1); 
				}
				else 
				{
					input_file = TFile::Open(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_w_error_data.root", i, run_lower, run_upper));
					if(!input_file) input_file = TFile::Open(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d.root", i, run_lower, run_upper));
					c = (TCanvas*)input_file->Get("canvas");
					if(!c) c = (TCanvas*)input_file->Get(Form("Layer%d_fakehitrate_w_error_and_trig_data", i));
					
				    c->Draw();
				    c->SaveAs(Form("FHR_from_run%d_to_run%d_L%d.png", run_lower, run_upper, i));
				    n_stave = number_staves[i];
				    vector<TGraph*> staves;
				    for(Int_t j=0; j < n_stave; j++) 
				    {
				    	TGraph* aux_stave = (TGraph*) c->GetPrimitive((Form("gr_L%d_stave%d_HS0", i, j)));
				    	staves.push_back(aux_stave);
				    }
	    			layer_name = Form("L%d", i);
	    			create_repository(i, n_stave, total_events, cumulative_label, label_event[i], staves, layer_name, final_comment, build_vcomments, build_layer_status, build_layer_collision);
	    			comments.push_back(build_vcomments);
	    			layer_status.push_back(build_layer_status);
	    			layer_collision.push_back(build_layer_collision);
	    			build_vcomments.clear();
	    			build_layer_status.clear();
	    			build_layer_collision.clear();
	    			cout << "repository operation on layer " << i << " concluded " << endl;
	    			staves.clear();			
				}
			}
			else
			{
				if(gSystem->AccessPathName(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_w_error_data_HSLower.root", i, run_lower, run_upper))&&gSystem->AccessPathName(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_HSLower.root", i, run_lower, run_upper)))
				{
					cout << "The file for layer " << i << " HS lower does not exist. Check the run extremes. " << endl;
					exit(1);
				}
				else
				{
					input_file = TFile::Open(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_w_error_data_HSLower.root", i, run_lower, run_upper));
					if(!input_file) input_file = TFile::Open(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_HSLower.root", i, run_lower, run_upper));
					c = (TCanvas*)input_file->Get("canvas_HSLower");
					if(!c) c = (TCanvas*)input_file->Get(Form("Layer%d_fakehitrate_w_error_and_trig_data_HSLower", i));
					
				    c->Draw();
				    c->SaveAs(Form("FHR_from_run%d_to_run%d_L%d_HSLower.png", run_lower, run_upper, i));
				    n_stave = number_staves[i];
				    vector<TGraph*> staves_lower;
				    for(Int_t j=0; j < n_stave; j++) 
				    {
					    TGraph* aux_stave = (TGraph*) c->GetPrimitive((Form("gr_L%d_stave%d_HS0", i, j)));
					    staves_lower.push_back(aux_stave);
				    }
				    layer_name = Form("L%dL", i);
				    create_repository(i, n_stave, total_events, cumulative_label, label_event[i], staves_lower, layer_name, final_comment, build_vcomments, build_layer_status, build_layer_collision);
			    	comments.push_back(build_vcomments);
			    	layer_status.push_back(build_layer_status);
			    	layer_collision.push_back(build_layer_collision);
			    	build_vcomments.clear();
			    	build_layer_status.clear();
			    	build_layer_collision.clear();
			    	cout << "repository operation on layer " << i << " lower concluded " << endl;
			    	staves_lower.clear();					
				}
				
				if(gSystem->AccessPathName(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_w_error_data_HSUpper.root", i, run_lower, run_upper))&&gSystem->AccessPathName(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_HSUpper.root", i, run_lower, run_upper)))
				{
				    cout << "The file for layer " << i << " HS upper does not exist. Check the run extremes. " << endl;
		    	    exit(1);
				}
				else
				{
					input_file = TFile::Open(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_w_error_data_HSUpper.root", i, run_lower, run_upper));
					if(!input_file) input_file = TFile::Open(Form("../Plots/Layer%d_fakehitrate_from_run%d_to_run%d_HSUpper.root", i, run_lower, run_upper));
					c = (TCanvas*)input_file->Get("canvas_HSUpper");
					if(!c) c = (TCanvas*)input_file->Get(Form("Layer%d_fakehitrate_w_error_and_trig_data_HSUpper", i));
					
				    c->Draw();
				    c->SaveAs(Form("FHR_from_run%d_to_run%d_L%d_HSUpper.png", run_lower, run_upper, i));
				    n_stave = number_staves[i];
				    vector<TGraph*> staves_upper;
				    for(Int_t j=0; j < n_stave; j++) 
				    {
				    	TGraph* aux_stave = (TGraph*) c->GetPrimitive((Form("gr_L%d_stave%d_HS1", i, j)));
				    	staves_upper.push_back(aux_stave);
				    }
				    layer_name = Form("L%dU", i);
				    create_repository(i, n_stave, total_events, cumulative_label, label_event[i], staves_upper, layer_name, final_comment, build_vcomments, build_layer_status, build_layer_collision);
				    comments.push_back(build_vcomments);
				    layer_status.push_back(build_layer_status);
				    layer_collision.push_back(build_layer_collision);
				    build_vcomments.clear();
				    build_layer_status.clear();
				    build_layer_collision.clear();
				    cout << "repository operation on layer " << i << " upper concluded " << endl;
				    staves_upper.clear();			
				}
			}
		}
		
		for(Int_t j = 0; j < total_events; j++)
		{
			for(Int_t i=0; i < 11; i++)
			{
				build_definitive_comments.append(comments[i][j]);
			}
			build_definitive_comments.pop_back();
			definitive_comments.push_back(build_definitive_comments);
			build_definitive_comments = "";
		}
		
		for(Int_t j=0; j < total_events; j++)
		{
			min_status = 10;
			Int_t layer_status_size = layer_status.size(); 
			for(Int_t i=0; i < layer_status_size; i++)
			{
				if(min_status > layer_status[i][j]) min_status = layer_status[i][j];
			}
			global_status.push_back(min_status);
		}
		Int_t global_status_size = global_status.size();
		for(Int_t i=0; i < global_status_size; i++)
		{
			if(global_status[i]==0) global_status_descriptive.push_back("expert check");
			else if(global_status[i]==1) global_status_descriptive.push_back("bad");
			else global_status_descriptive.push_back("good");
		}
	
		for(Int_t i=0; i < 11; i++)
		{
			for(Int_t j=0; j < total_events; j++)
			{
				if(layer_status[i][j]==0) build_layer_status_descriptive.push_back("expert check");
				else if(layer_status[i][j]==1) build_layer_status_descriptive.push_back("bad");
				else if(layer_status[i][j]==2) build_layer_status_descriptive.push_back("good");
				else if(layer_status[i][j]==3) build_layer_status_descriptive.push_back("N/A");
			}
			layer_status_descriptive.push_back(build_layer_status_descriptive);
			build_layer_status_descriptive.clear();
		}
	
		for(Int_t j=0; j < total_events; j++)
		{
			Int_t counter = 0;
			collision_type = layer_collision[0][j]; 
			Int_t layer_collision_size = layer_collision.size();
			for(Int_t i=0; i < layer_collision_size; i++)
			{
				if(collision_type!=layer_collision[i][j]) counter++;
			}
			if(counter==0) global_collision.push_back(layer_collision[0][j]); 
			else global_collision.push_back(0);
		}
		Int_t global_collision_size = global_collision.size();
		for(Int_t i=0; i < global_collision_size; i++)
		{
			if (global_collision[i]==0) global_collision_descriptive.push_back("expert check");
			else if(global_collision[i]==1) global_collision_descriptive.push_back("cosmic");
			else if(global_collision[i]==2) global_collision_descriptive.push_back("collisions pp");
			else if(global_collision[i]==3) global_collision_descriptive.push_back("collisions Pb-Pb");
			else if(global_collision[i]==4) global_collision_descriptive.push_back("emulated pp");
			else if(global_collision[i]==5) global_collision_descriptive.push_back("emulated Pb-Pb");
		}
	
		ofstream output;
		output.open(Form("repository_from_run%d_to_run%d.txt", run_lower, run_upper));
		if(output.fail())
		{
			cout << "Error in repository file creation " << endl;
			exit(1);
		}
		
		for(Int_t i=0; i < total_events; i++) 
		{
			output << cumulative_label[i] << "\t" << global_collision_descriptive[i] << "\t" << global_status_descriptive[i] << "\t";
			for(Int_t j=0; j < 11; j++)
			{
				output << layer_status_descriptive[j][i] << "\t";
			}
			output << definitive_comments[i] << endl;
		}
	}
	
}
