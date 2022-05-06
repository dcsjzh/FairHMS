#ifndef __HITTINGSET_H__
#define __HITTINGSET_H__

#include <map>
#include <set>
#include <vector>

using namespace std;

class HSApprox {
  public:
    //S is the collection of sets, specified as vector of sets
   //HS. static void get_hs_approximation2(const vector<set<size_t> >& S, vector< size_t >& HS);
    template < typename T >
    static void get_hs_approximation(const vector<set< T> >& S, vector< T >& HS);        
};

template <typename T>
void HSApprox::get_hs_approximation(const vector< set < T > >& S, vector < T >& HS)
{
  map < T, set < size_t >  > cmap;
  for(size_t i = 0; i < S.size() ; i++)
  {
    class set< T > :: iterator it;
    for(it = S[i].begin(); it != S[i].end(); ++it)
    {
      T el = (*it);
      class map < T, set < size_t > > :: iterator Xit = cmap.find(el);

      //if this element has been seen before
      if(Xit != cmap.end())
        (Xit->second).insert(i);
      else //creat a new entry for this element
      {
	set < size_t > NSet;
	NSet.insert(i);
        cmap[el] = NSet;
      }
    }
  }

  //now output the hitting set
  set < size_t > uncovered_sets;
  for(size_t i = 0; i < S.size(); i++)
    uncovered_sets.insert(i);

  while(!uncovered_sets.empty())
  {
  //  if(HS.size()==5){
    //    return;}
    //pick the element that hits maximum number of sets
    class map < T, set < size_t > > :: iterator Xit;
    int smax = -1;
    T emax;
    set < size_t > Hset;
    for(Xit = cmap.begin(); Xit != cmap.end(); ++Xit)
    {
      T el = Xit->first;
      int hits = 0;


      set < size_t > :: iterator itA;

      //do the right count - count only uncovered sets
      for(itA = Xit->second.begin(); itA != Xit->second.end(); itA++)
      {
        if(uncovered_sets.find(*itA) != uncovered_sets.end())
	  hits++;
      }

      if(hits > smax)
      {
        smax = hits;
        emax = el;
        Hset = Xit->second;
      }
    }

    assert(smax > 0);

    set < size_t > :: iterator itA;
    //now delete all uncovered sets hit by element
    //add el to the hitting set and delete from cmap
    for(itA = Hset.begin(); itA !=  Hset.end(); itA++)
    {
      set < size_t > :: iterator it = uncovered_sets.find(*itA);
      
      if(it != uncovered_sets.end()) 
        uncovered_sets.erase(it);
    }

    HS.push_back(emax);
  }
}


/*
void HSApprox::get_hs_approximation2(const vector< set <size_t> >& S, vector < size_t >& HS)
{
  map < size_t, set < size_t >  > cmap;
  for(size_t i = 0; i < S.size() ; i++)
  {
    class set< size_t > :: iterator it;
    for(it = S[i].begin(); it != S[i].end(); ++it)
    {
      size_t el = (*it);
      class map < size_t, set < size_t > > :: iterator Xit = cmap.find(el);

      //if this element has been seen before
      if(Xit != cmap.end())
        (Xit->second).insert(i);
      else //creat a new entry for this element
      {
	set < size_t > NSet;
	NSet.insert(i);
        cmap[el] = NSet;
      }
    }
  }

  //now output the hitting set
  set < size_t > uncovered_sets;
  for(size_t i = 0; i < S.size(); i++)
    uncovered_sets.insert(i);

  while(!uncovered_sets.empty())
  {
    if(HS.size()==20){
        return;
    }
    //pick the element that hits maximum number of sets
    class map < size_t, set < size_t > > :: iterator Xit;
    int smax = -1;
    size_t emax;
    set < size_t > Hset;
    for(Xit = cmap.begin(); Xit != cmap.end(); ++Xit)
    {
      size_t el = Xit->first;
      int hits = 0;


      set < size_t > :: iterator itA;

      //do the right count - count only uncovered sets
      for(itA = Xit->second.begin(); itA != Xit->second.end(); itA++)
      {
        if(uncovered_sets.find(*itA) != uncovered_sets.end())
	  hits++;
      }

      if(hits > smax)
      {
        smax = hits;
        emax = el;
        Hset = Xit->second;
      }
    }

    assert(smax > 0);

    set < size_t > :: iterator itA;
    //now delete all uncovered sets hit by element
    //add el to the hitting set and delete from cmap
    for(itA = Hset.begin(); itA !=  Hset.end(); itA++)
    {
      set < size_t > :: iterator it = uncovered_sets.find(*itA);
      
      if(it != uncovered_sets.end()) 
        uncovered_sets.erase(it);
    }

    HS.push_back(emax);
  }
}
*/
#endif
