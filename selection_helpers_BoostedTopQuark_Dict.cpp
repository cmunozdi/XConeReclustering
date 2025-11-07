// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME selection_helpers_BoostedTopQuark_Dict
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "ROOT/RConfig.hxx"
#include "TClass.h"
#include "TDictAttributeMap.h"
#include "TInterpreter.h"
#include "TROOT.h"
#include "TBuffer.h"
#include "TMemberInspector.h"
#include "TInterpreter.h"
#include "TVirtualMutex.h"
#include "TError.h"

#ifndef G__ROOT
#define G__ROOT
#endif

#include "RtypesImp.h"
#include "TIsAProxy.h"
#include "TFileMergeInfo.h"
#include <algorithm>
#include "TCollectionProxyInfo.h"
/*******************************************************************/

#include "TDataMember.h"

// Header files passed as explicit arguments
#include "selection_helpers_BoostedTopQuark.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *JetAntikTReclus_Dictionary();
   static void JetAntikTReclus_TClassManip(TClass*);
   static void *new_JetAntikTReclus(void *p = nullptr);
   static void *newArray_JetAntikTReclus(Long_t size, void *p);
   static void delete_JetAntikTReclus(void *p);
   static void deleteArray_JetAntikTReclus(void *p);
   static void destruct_JetAntikTReclus(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::JetAntikTReclus*)
   {
      ::JetAntikTReclus *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(::JetAntikTReclus));
      static ::ROOT::TGenericClassInfo 
         instance("JetAntikTReclus", "selection_helpers_BoostedTopQuark.h", 194,
                  typeid(::JetAntikTReclus), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &JetAntikTReclus_Dictionary, isa_proxy, 4,
                  sizeof(::JetAntikTReclus) );
      instance.SetNew(&new_JetAntikTReclus);
      instance.SetNewArray(&newArray_JetAntikTReclus);
      instance.SetDelete(&delete_JetAntikTReclus);
      instance.SetDeleteArray(&deleteArray_JetAntikTReclus);
      instance.SetDestructor(&destruct_JetAntikTReclus);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::JetAntikTReclus*)
   {
      return GenerateInitInstanceLocal(static_cast<::JetAntikTReclus*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::JetAntikTReclus*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *JetAntikTReclus_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const ::JetAntikTReclus*>(nullptr))->GetClass();
      JetAntikTReclus_TClassManip(theClass);
   return theClass;
   }

   static void JetAntikTReclus_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   static void *new_Lepton(void *p = nullptr);
   static void *newArray_Lepton(Long_t size, void *p);
   static void delete_Lepton(void *p);
   static void deleteArray_Lepton(void *p);
   static void destruct_Lepton(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::Lepton*)
   {
      ::Lepton *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::Lepton >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("Lepton", ::Lepton::Class_Version(), "selection_helpers_BoostedTopQuark.h", 240,
                  typeid(::Lepton), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::Lepton::Dictionary, isa_proxy, 4,
                  sizeof(::Lepton) );
      instance.SetNew(&new_Lepton);
      instance.SetNewArray(&newArray_Lepton);
      instance.SetDelete(&delete_Lepton);
      instance.SetDeleteArray(&deleteArray_Lepton);
      instance.SetDestructor(&destruct_Lepton);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::Lepton*)
   {
      return GenerateInitInstanceLocal(static_cast<::Lepton*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::Lepton*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_JetReclus(void *p = nullptr);
   static void *newArray_JetReclus(Long_t size, void *p);
   static void delete_JetReclus(void *p);
   static void deleteArray_JetReclus(void *p);
   static void destruct_JetReclus(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::JetReclus*)
   {
      ::JetReclus *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::JetReclus >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("JetReclus", ::JetReclus::Class_Version(), "selection_helpers_BoostedTopQuark.h", 728,
                  typeid(::JetReclus), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::JetReclus::Dictionary, isa_proxy, 4,
                  sizeof(::JetReclus) );
      instance.SetNew(&new_JetReclus);
      instance.SetNewArray(&newArray_JetReclus);
      instance.SetDelete(&delete_JetReclus);
      instance.SetDeleteArray(&deleteArray_JetReclus);
      instance.SetDestructor(&destruct_JetReclus);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::JetReclus*)
   {
      return GenerateInitInstanceLocal(static_cast<::JetReclus*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::JetReclus*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_TopJetReclus(void *p = nullptr);
   static void *newArray_TopJetReclus(Long_t size, void *p);
   static void delete_TopJetReclus(void *p);
   static void deleteArray_TopJetReclus(void *p);
   static void destruct_TopJetReclus(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::TopJetReclus*)
   {
      ::TopJetReclus *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::TopJetReclus >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("TopJetReclus", ::TopJetReclus::Class_Version(), "selection_helpers_BoostedTopQuark.h", 749,
                  typeid(::TopJetReclus), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::TopJetReclus::Dictionary, isa_proxy, 4,
                  sizeof(::TopJetReclus) );
      instance.SetNew(&new_TopJetReclus);
      instance.SetNewArray(&newArray_TopJetReclus);
      instance.SetDelete(&delete_TopJetReclus);
      instance.SetDeleteArray(&deleteArray_TopJetReclus);
      instance.SetDestructor(&destruct_TopJetReclus);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::TopJetReclus*)
   {
      return GenerateInitInstanceLocal(static_cast<::TopJetReclus*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::TopJetReclus*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_XConeReclusteredJets(void *p = nullptr);
   static void *newArray_XConeReclusteredJets(Long_t size, void *p);
   static void delete_XConeReclusteredJets(void *p);
   static void deleteArray_XConeReclusteredJets(void *p);
   static void destruct_XConeReclusteredJets(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::XConeReclusteredJets*)
   {
      ::XConeReclusteredJets *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::XConeReclusteredJets >(nullptr);
      static ::ROOT::TGenericClassInfo 
         instance("XConeReclusteredJets", ::XConeReclusteredJets::Class_Version(), "selection_helpers_BoostedTopQuark.h", 770,
                  typeid(::XConeReclusteredJets), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::XConeReclusteredJets::Dictionary, isa_proxy, 4,
                  sizeof(::XConeReclusteredJets) );
      instance.SetNew(&new_XConeReclusteredJets);
      instance.SetNewArray(&newArray_XConeReclusteredJets);
      instance.SetDelete(&delete_XConeReclusteredJets);
      instance.SetDeleteArray(&deleteArray_XConeReclusteredJets);
      instance.SetDestructor(&destruct_XConeReclusteredJets);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::XConeReclusteredJets*)
   {
      return GenerateInitInstanceLocal(static_cast<::XConeReclusteredJets*>(nullptr));
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const ::XConeReclusteredJets*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

//______________________________________________________________________________
atomic_TClass_ptr Lepton::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *Lepton::Class_Name()
{
   return "Lepton";
}

//______________________________________________________________________________
const char *Lepton::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Lepton*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int Lepton::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::Lepton*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Lepton::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Lepton*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Lepton::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::Lepton*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr JetReclus::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *JetReclus::Class_Name()
{
   return "JetReclus";
}

//______________________________________________________________________________
const char *JetReclus::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::JetReclus*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int JetReclus::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::JetReclus*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *JetReclus::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::JetReclus*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *JetReclus::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::JetReclus*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr TopJetReclus::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *TopJetReclus::Class_Name()
{
   return "TopJetReclus";
}

//______________________________________________________________________________
const char *TopJetReclus::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TopJetReclus*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int TopJetReclus::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::TopJetReclus*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *TopJetReclus::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TopJetReclus*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *TopJetReclus::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::TopJetReclus*)nullptr)->GetClass(); }
   return fgIsA;
}

//______________________________________________________________________________
atomic_TClass_ptr XConeReclusteredJets::fgIsA(nullptr);  // static to hold class pointer

//______________________________________________________________________________
const char *XConeReclusteredJets::Class_Name()
{
   return "XConeReclusteredJets";
}

//______________________________________________________________________________
const char *XConeReclusteredJets::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::XConeReclusteredJets*)nullptr)->GetImplFileName();
}

//______________________________________________________________________________
int XConeReclusteredJets::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::XConeReclusteredJets*)nullptr)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *XConeReclusteredJets::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::XConeReclusteredJets*)nullptr)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *XConeReclusteredJets::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::XConeReclusteredJets*)nullptr)->GetClass(); }
   return fgIsA;
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_JetAntikTReclus(void *p) {
      return  p ? new(p) ::JetAntikTReclus : new ::JetAntikTReclus;
   }
   static void *newArray_JetAntikTReclus(Long_t nElements, void *p) {
      return p ? new(p) ::JetAntikTReclus[nElements] : new ::JetAntikTReclus[nElements];
   }
   // Wrapper around operator delete
   static void delete_JetAntikTReclus(void *p) {
      delete (static_cast<::JetAntikTReclus*>(p));
   }
   static void deleteArray_JetAntikTReclus(void *p) {
      delete [] (static_cast<::JetAntikTReclus*>(p));
   }
   static void destruct_JetAntikTReclus(void *p) {
      typedef ::JetAntikTReclus current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::JetAntikTReclus

//______________________________________________________________________________
void Lepton::Streamer(TBuffer &R__b)
{
   // Stream an object of class Lepton.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(Lepton::Class(),this);
   } else {
      R__b.WriteClassBuffer(Lepton::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_Lepton(void *p) {
      return  p ? new(p) ::Lepton : new ::Lepton;
   }
   static void *newArray_Lepton(Long_t nElements, void *p) {
      return p ? new(p) ::Lepton[nElements] : new ::Lepton[nElements];
   }
   // Wrapper around operator delete
   static void delete_Lepton(void *p) {
      delete (static_cast<::Lepton*>(p));
   }
   static void deleteArray_Lepton(void *p) {
      delete [] (static_cast<::Lepton*>(p));
   }
   static void destruct_Lepton(void *p) {
      typedef ::Lepton current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::Lepton

//______________________________________________________________________________
void JetReclus::Streamer(TBuffer &R__b)
{
   // Stream an object of class JetReclus.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(JetReclus::Class(),this);
   } else {
      R__b.WriteClassBuffer(JetReclus::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_JetReclus(void *p) {
      return  p ? new(p) ::JetReclus : new ::JetReclus;
   }
   static void *newArray_JetReclus(Long_t nElements, void *p) {
      return p ? new(p) ::JetReclus[nElements] : new ::JetReclus[nElements];
   }
   // Wrapper around operator delete
   static void delete_JetReclus(void *p) {
      delete (static_cast<::JetReclus*>(p));
   }
   static void deleteArray_JetReclus(void *p) {
      delete [] (static_cast<::JetReclus*>(p));
   }
   static void destruct_JetReclus(void *p) {
      typedef ::JetReclus current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::JetReclus

//______________________________________________________________________________
void TopJetReclus::Streamer(TBuffer &R__b)
{
   // Stream an object of class TopJetReclus.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(TopJetReclus::Class(),this);
   } else {
      R__b.WriteClassBuffer(TopJetReclus::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_TopJetReclus(void *p) {
      return  p ? new(p) ::TopJetReclus : new ::TopJetReclus;
   }
   static void *newArray_TopJetReclus(Long_t nElements, void *p) {
      return p ? new(p) ::TopJetReclus[nElements] : new ::TopJetReclus[nElements];
   }
   // Wrapper around operator delete
   static void delete_TopJetReclus(void *p) {
      delete (static_cast<::TopJetReclus*>(p));
   }
   static void deleteArray_TopJetReclus(void *p) {
      delete [] (static_cast<::TopJetReclus*>(p));
   }
   static void destruct_TopJetReclus(void *p) {
      typedef ::TopJetReclus current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::TopJetReclus

//______________________________________________________________________________
void XConeReclusteredJets::Streamer(TBuffer &R__b)
{
   // Stream an object of class XConeReclusteredJets.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(XConeReclusteredJets::Class(),this);
   } else {
      R__b.WriteClassBuffer(XConeReclusteredJets::Class(),this);
   }
}

namespace ROOT {
   // Wrappers around operator new
   static void *new_XConeReclusteredJets(void *p) {
      return  p ? new(p) ::XConeReclusteredJets : new ::XConeReclusteredJets;
   }
   static void *newArray_XConeReclusteredJets(Long_t nElements, void *p) {
      return p ? new(p) ::XConeReclusteredJets[nElements] : new ::XConeReclusteredJets[nElements];
   }
   // Wrapper around operator delete
   static void delete_XConeReclusteredJets(void *p) {
      delete (static_cast<::XConeReclusteredJets*>(p));
   }
   static void deleteArray_XConeReclusteredJets(void *p) {
      delete [] (static_cast<::XConeReclusteredJets*>(p));
   }
   static void destruct_XConeReclusteredJets(void *p) {
      typedef ::XConeReclusteredJets current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class ::XConeReclusteredJets

namespace ROOT {
   static TClass *vectorlEintgR_Dictionary();
   static void vectorlEintgR_TClassManip(TClass*);
   static void *new_vectorlEintgR(void *p = nullptr);
   static void *newArray_vectorlEintgR(Long_t size, void *p);
   static void delete_vectorlEintgR(void *p);
   static void deleteArray_vectorlEintgR(void *p);
   static void destruct_vectorlEintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<int>*)
   {
      vector<int> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<int>", -2, "vector", 423,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<int>","std::vector<int, std::allocator<int> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<int>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<int>*>(nullptr))->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete (static_cast<vector<int>*>(p));
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] (static_cast<vector<int>*>(p));
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<int>

namespace ROOT {
   static TClass *vectorlEfloatgR_Dictionary();
   static void vectorlEfloatgR_TClassManip(TClass*);
   static void *new_vectorlEfloatgR(void *p = nullptr);
   static void *newArray_vectorlEfloatgR(Long_t size, void *p);
   static void delete_vectorlEfloatgR(void *p);
   static void deleteArray_vectorlEfloatgR(void *p);
   static void destruct_vectorlEfloatgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<float>*)
   {
      vector<float> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<float>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<float>", -2, "vector", 423,
                  typeid(vector<float>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEfloatgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<float>) );
      instance.SetNew(&new_vectorlEfloatgR);
      instance.SetNewArray(&newArray_vectorlEfloatgR);
      instance.SetDelete(&delete_vectorlEfloatgR);
      instance.SetDeleteArray(&deleteArray_vectorlEfloatgR);
      instance.SetDestructor(&destruct_vectorlEfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<float> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<float>","std::vector<float, std::allocator<float> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<float>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<float>*>(nullptr))->GetClass();
      vectorlEfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEfloatgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<float> : new vector<float>;
   }
   static void *newArray_vectorlEfloatgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<float>[nElements] : new vector<float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEfloatgR(void *p) {
      delete (static_cast<vector<float>*>(p));
   }
   static void deleteArray_vectorlEfloatgR(void *p) {
      delete [] (static_cast<vector<float>*>(p));
   }
   static void destruct_vectorlEfloatgR(void *p) {
      typedef vector<float> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<float>

namespace ROOT {
   static TClass *vectorlEboolgR_Dictionary();
   static void vectorlEboolgR_TClassManip(TClass*);
   static void *new_vectorlEboolgR(void *p = nullptr);
   static void *newArray_vectorlEboolgR(Long_t size, void *p);
   static void delete_vectorlEboolgR(void *p);
   static void deleteArray_vectorlEboolgR(void *p);
   static void destruct_vectorlEboolgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<bool>*)
   {
      vector<bool> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<bool>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<bool>", -2, "vector", 690,
                  typeid(vector<bool>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEboolgR_Dictionary, isa_proxy, 0,
                  sizeof(vector<bool>) );
      instance.SetNew(&new_vectorlEboolgR);
      instance.SetNewArray(&newArray_vectorlEboolgR);
      instance.SetDelete(&delete_vectorlEboolgR);
      instance.SetDeleteArray(&deleteArray_vectorlEboolgR);
      instance.SetDestructor(&destruct_vectorlEboolgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<bool> >()));

      instance.AdoptAlternate(::ROOT::AddClassAlternate("vector<bool>","std::vector<bool, std::allocator<bool> >"));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal(static_cast<const vector<bool>*>(nullptr)); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEboolgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal(static_cast<const vector<bool>*>(nullptr))->GetClass();
      vectorlEboolgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEboolgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEboolgR(void *p) {
      return  p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<bool> : new vector<bool>;
   }
   static void *newArray_vectorlEboolgR(Long_t nElements, void *p) {
      return p ? ::new(static_cast<::ROOT::Internal::TOperatorNewHelper*>(p)) vector<bool>[nElements] : new vector<bool>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEboolgR(void *p) {
      delete (static_cast<vector<bool>*>(p));
   }
   static void deleteArray_vectorlEboolgR(void *p) {
      delete [] (static_cast<vector<bool>*>(p));
   }
   static void destruct_vectorlEboolgR(void *p) {
      typedef vector<bool> current_t;
      (static_cast<current_t*>(p))->~current_t();
   }
} // end of namespace ROOT for class vector<bool>

namespace {
  void TriggerDictionaryInitialization_selection_helpers_BoostedTopQuark_Dict_Impl() {
    static const char* headers[] = {
"selection_helpers_BoostedTopQuark.h",
nullptr
    };
    static const char* includePaths[] = {
"/cvmfs/cms.cern.ch/el9_amd64_gcc12/lcg/root/6.32.11-a1d23f0c91e1aeff81bd016882a5fb71/include/",
"/afs/cern.ch/user/c/cmunozdi/CMSSW_15_0_3_patch1/src/XConeReclustering/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "selection_helpers_BoostedTopQuark_Dict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
struct __attribute__((annotate("$clingAutoload$selection_helpers_BoostedTopQuark.h")))  JetAntikTReclus;
struct __attribute__((annotate("$clingAutoload$selection_helpers_BoostedTopQuark.h")))  Lepton;
struct __attribute__((annotate("$clingAutoload$selection_helpers_BoostedTopQuark.h")))  JetReclus;
struct __attribute__((annotate("$clingAutoload$selection_helpers_BoostedTopQuark.h")))  TopJetReclus;
struct __attribute__((annotate("$clingAutoload$selection_helpers_BoostedTopQuark.h")))  XConeReclusteredJets;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "selection_helpers_BoostedTopQuark_Dict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "selection_helpers_BoostedTopQuark.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"JetAntikTReclus", payloadCode, "@",
"JetReclus", payloadCode, "@",
"Lepton", payloadCode, "@",
"TopJetReclus", payloadCode, "@",
"XConeReclusteredJets", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("selection_helpers_BoostedTopQuark_Dict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_selection_helpers_BoostedTopQuark_Dict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_selection_helpers_BoostedTopQuark_Dict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_selection_helpers_BoostedTopQuark_Dict() {
  TriggerDictionaryInitialization_selection_helpers_BoostedTopQuark_Dict_Impl();
}
