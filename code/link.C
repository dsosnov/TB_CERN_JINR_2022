// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME link
#define R__NO_DEPRECATION

/*******************************************************************/
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#define G__DICTIONARY
#include "RConfig.h"
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
#include "link.h"

// Header files passed via #pragma extra_include

// The generated code does not explicitly qualify STL entities
namespace std {} using namespace std;

namespace ROOT {
   static TClass *vectorlEvectorlEunsignedsPintgRsPgR_Dictionary();
   static void vectorlEvectorlEunsignedsPintgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEunsignedsPintgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEunsignedsPintgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEunsignedsPintgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEunsignedsPintgRsPgR(void *p);
   static void destruct_vectorlEvectorlEunsignedsPintgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<unsigned int> >*)
   {
      vector<vector<unsigned int> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<unsigned int> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<unsigned int> >", -2, "c++/v1/vector", 471,
                  typeid(vector<vector<unsigned int> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEunsignedsPintgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<unsigned int> >) );
      instance.SetNew(&new_vectorlEvectorlEunsignedsPintgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEunsignedsPintgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEunsignedsPintgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEunsignedsPintgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEunsignedsPintgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<unsigned int> > >()));

      ::ROOT::AddClassAlternate("vector<vector<unsigned int> >","std::__1::vector<std::__1::vector<unsigned int, std::__1::allocator<unsigned int> >, std::__1::allocator<std::__1::vector<unsigned int, std::__1::allocator<unsigned int> > > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vector<unsigned int> >*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEunsignedsPintgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<unsigned int> >*)nullptr)->GetClass();
      vectorlEvectorlEunsignedsPintgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEunsignedsPintgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEunsignedsPintgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<unsigned int> > : new vector<vector<unsigned int> >;
   }
   static void *newArray_vectorlEvectorlEunsignedsPintgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<unsigned int> >[nElements] : new vector<vector<unsigned int> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEunsignedsPintgRsPgR(void *p) {
      delete ((vector<vector<unsigned int> >*)p);
   }
   static void deleteArray_vectorlEvectorlEunsignedsPintgRsPgR(void *p) {
      delete [] ((vector<vector<unsigned int> >*)p);
   }
   static void destruct_vectorlEvectorlEunsignedsPintgRsPgR(void *p) {
      typedef vector<vector<unsigned int> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<unsigned int> >

namespace ROOT {
   static TClass *vectorlEvectorlEstringgRsPgR_Dictionary();
   static void vectorlEvectorlEstringgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEstringgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEstringgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEstringgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEstringgRsPgR(void *p);
   static void destruct_vectorlEvectorlEstringgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<string> >*)
   {
      vector<vector<string> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<string> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<string> >", -2, "c++/v1/vector", 471,
                  typeid(vector<vector<string> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEstringgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<string> >) );
      instance.SetNew(&new_vectorlEvectorlEstringgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEstringgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEstringgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEstringgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEstringgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<string> > >()));

      ::ROOT::AddClassAlternate("vector<vector<string> >","std::__1::vector<std::__1::vector<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char> >, std::__1::allocator<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char> > > >, std::__1::allocator<std::__1::vector<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char> >, std::__1::allocator<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char> > > > > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vector<string> >*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEstringgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<string> >*)nullptr)->GetClass();
      vectorlEvectorlEstringgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEstringgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEstringgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<string> > : new vector<vector<string> >;
   }
   static void *newArray_vectorlEvectorlEstringgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<string> >[nElements] : new vector<vector<string> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEstringgRsPgR(void *p) {
      delete ((vector<vector<string> >*)p);
   }
   static void deleteArray_vectorlEvectorlEstringgRsPgR(void *p) {
      delete [] ((vector<vector<string> >*)p);
   }
   static void destruct_vectorlEvectorlEstringgRsPgR(void *p) {
      typedef vector<vector<string> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<string> >

namespace ROOT {
   static TClass *vectorlEvectorlEshortgRsPgR_Dictionary();
   static void vectorlEvectorlEshortgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEshortgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEshortgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEshortgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEshortgRsPgR(void *p);
   static void destruct_vectorlEvectorlEshortgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<short> >*)
   {
      vector<vector<short> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<short> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<short> >", -2, "c++/v1/vector", 471,
                  typeid(vector<vector<short> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEshortgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<short> >) );
      instance.SetNew(&new_vectorlEvectorlEshortgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEshortgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEshortgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEshortgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEshortgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<short> > >()));

      ::ROOT::AddClassAlternate("vector<vector<short> >","std::__1::vector<std::__1::vector<short, std::__1::allocator<short> >, std::__1::allocator<std::__1::vector<short, std::__1::allocator<short> > > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vector<short> >*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEshortgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<short> >*)nullptr)->GetClass();
      vectorlEvectorlEshortgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEshortgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEshortgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<short> > : new vector<vector<short> >;
   }
   static void *newArray_vectorlEvectorlEshortgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<short> >[nElements] : new vector<vector<short> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEshortgRsPgR(void *p) {
      delete ((vector<vector<short> >*)p);
   }
   static void deleteArray_vectorlEvectorlEshortgRsPgR(void *p) {
      delete [] ((vector<vector<short> >*)p);
   }
   static void destruct_vectorlEvectorlEshortgRsPgR(void *p) {
      typedef vector<vector<short> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<short> >

namespace ROOT {
   static TClass *vectorlEvectorlEintgRsPgR_Dictionary();
   static void vectorlEvectorlEintgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEintgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEintgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEintgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEintgRsPgR(void *p);
   static void destruct_vectorlEvectorlEintgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<int> >*)
   {
      vector<vector<int> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<int> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<int> >", -2, "c++/v1/vector", 471,
                  typeid(vector<vector<int> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEintgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<int> >) );
      instance.SetNew(&new_vectorlEvectorlEintgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEintgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEintgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEintgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEintgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<int> > >()));

      ::ROOT::AddClassAlternate("vector<vector<int> >","std::__1::vector<std::__1::vector<int, std::__1::allocator<int> >, std::__1::allocator<std::__1::vector<int, std::__1::allocator<int> > > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vector<int> >*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEintgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<int> >*)nullptr)->GetClass();
      vectorlEvectorlEintgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEintgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEintgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<int> > : new vector<vector<int> >;
   }
   static void *newArray_vectorlEvectorlEintgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<int> >[nElements] : new vector<vector<int> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEintgRsPgR(void *p) {
      delete ((vector<vector<int> >*)p);
   }
   static void deleteArray_vectorlEvectorlEintgRsPgR(void *p) {
      delete [] ((vector<vector<int> >*)p);
   }
   static void destruct_vectorlEvectorlEintgRsPgR(void *p) {
      typedef vector<vector<int> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<int> >

namespace ROOT {
   static TClass *vectorlEvectorlEfloatgRsPgR_Dictionary();
   static void vectorlEvectorlEfloatgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEfloatgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEfloatgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEfloatgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEfloatgRsPgR(void *p);
   static void destruct_vectorlEvectorlEfloatgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<float> >*)
   {
      vector<vector<float> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<float> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<float> >", -2, "c++/v1/vector", 471,
                  typeid(vector<vector<float> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEfloatgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<float> >) );
      instance.SetNew(&new_vectorlEvectorlEfloatgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEfloatgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEfloatgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEfloatgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEfloatgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<float> > >()));

      ::ROOT::AddClassAlternate("vector<vector<float> >","std::__1::vector<std::__1::vector<float, std::__1::allocator<float> >, std::__1::allocator<std::__1::vector<float, std::__1::allocator<float> > > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vector<float> >*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEfloatgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<float> >*)nullptr)->GetClass();
      vectorlEvectorlEfloatgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEfloatgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEfloatgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<float> > : new vector<vector<float> >;
   }
   static void *newArray_vectorlEvectorlEfloatgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<float> >[nElements] : new vector<vector<float> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEfloatgRsPgR(void *p) {
      delete ((vector<vector<float> >*)p);
   }
   static void deleteArray_vectorlEvectorlEfloatgRsPgR(void *p) {
      delete [] ((vector<vector<float> >*)p);
   }
   static void destruct_vectorlEvectorlEfloatgRsPgR(void *p) {
      typedef vector<vector<float> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<float> >

namespace ROOT {
   static TClass *vectorlEvectorlEdoublegRsPgR_Dictionary();
   static void vectorlEvectorlEdoublegRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEdoublegRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEdoublegRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEdoublegRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEdoublegRsPgR(void *p);
   static void destruct_vectorlEvectorlEdoublegRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<double> >*)
   {
      vector<vector<double> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<double> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<double> >", -2, "c++/v1/vector", 471,
                  typeid(vector<vector<double> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEdoublegRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<double> >) );
      instance.SetNew(&new_vectorlEvectorlEdoublegRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEdoublegRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEdoublegRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEdoublegRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEdoublegRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<double> > >()));

      ::ROOT::AddClassAlternate("vector<vector<double> >","std::__1::vector<std::__1::vector<double, std::__1::allocator<double> >, std::__1::allocator<std::__1::vector<double, std::__1::allocator<double> > > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vector<double> >*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEdoublegRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<double> >*)nullptr)->GetClass();
      vectorlEvectorlEdoublegRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEdoublegRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEdoublegRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<double> > : new vector<vector<double> >;
   }
   static void *newArray_vectorlEvectorlEdoublegRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<double> >[nElements] : new vector<vector<double> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEdoublegRsPgR(void *p) {
      delete ((vector<vector<double> >*)p);
   }
   static void deleteArray_vectorlEvectorlEdoublegRsPgR(void *p) {
      delete [] ((vector<vector<double> >*)p);
   }
   static void destruct_vectorlEvectorlEdoublegRsPgR(void *p) {
      typedef vector<vector<double> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<double> >

namespace ROOT {
   static TClass *vectorlEvectorlEchargRsPgR_Dictionary();
   static void vectorlEvectorlEchargRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEchargRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEchargRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEchargRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEchargRsPgR(void *p);
   static void destruct_vectorlEvectorlEchargRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<char> >*)
   {
      vector<vector<char> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<char> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<char> >", -2, "c++/v1/vector", 471,
                  typeid(vector<vector<char> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEchargRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<char> >) );
      instance.SetNew(&new_vectorlEvectorlEchargRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEchargRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEchargRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEchargRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEchargRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<char> > >()));

      ::ROOT::AddClassAlternate("vector<vector<char> >","std::__1::vector<std::__1::vector<char, std::__1::allocator<char> >, std::__1::allocator<std::__1::vector<char, std::__1::allocator<char> > > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vector<char> >*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEchargRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<char> >*)nullptr)->GetClass();
      vectorlEvectorlEchargRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEchargRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEchargRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<char> > : new vector<vector<char> >;
   }
   static void *newArray_vectorlEvectorlEchargRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<char> >[nElements] : new vector<vector<char> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEchargRsPgR(void *p) {
      delete ((vector<vector<char> >*)p);
   }
   static void deleteArray_vectorlEvectorlEchargRsPgR(void *p) {
      delete [] ((vector<vector<char> >*)p);
   }
   static void destruct_vectorlEvectorlEchargRsPgR(void *p) {
      typedef vector<vector<char> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<char> >

namespace ROOT {
   static TClass *vectorlEvectorlEapvTrackgRsPgR_Dictionary();
   static void vectorlEvectorlEapvTrackgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEapvTrackgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEapvTrackgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEapvTrackgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEapvTrackgRsPgR(void *p);
   static void destruct_vectorlEvectorlEapvTrackgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<apvTrack> >*)
   {
      vector<vector<apvTrack> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<apvTrack> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<apvTrack> >", -2, "c++/v1/vector", 471,
                  typeid(vector<vector<apvTrack> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEapvTrackgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<apvTrack> >) );
      instance.SetNew(&new_vectorlEvectorlEapvTrackgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEapvTrackgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEapvTrackgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEapvTrackgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEapvTrackgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<apvTrack> > >()));

      ::ROOT::AddClassAlternate("vector<vector<apvTrack> >","std::__1::vector<std::__1::vector<apvTrack, std::__1::allocator<apvTrack> >, std::__1::allocator<std::__1::vector<apvTrack, std::__1::allocator<apvTrack> > > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vector<apvTrack> >*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEapvTrackgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<apvTrack> >*)nullptr)->GetClass();
      vectorlEvectorlEapvTrackgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEapvTrackgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEapvTrackgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<apvTrack> > : new vector<vector<apvTrack> >;
   }
   static void *newArray_vectorlEvectorlEapvTrackgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<apvTrack> >[nElements] : new vector<vector<apvTrack> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEapvTrackgRsPgR(void *p) {
      delete ((vector<vector<apvTrack> >*)p);
   }
   static void deleteArray_vectorlEvectorlEapvTrackgRsPgR(void *p) {
      delete [] ((vector<vector<apvTrack> >*)p);
   }
   static void destruct_vectorlEvectorlEapvTrackgRsPgR(void *p) {
      typedef vector<vector<apvTrack> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<apvTrack> >

namespace ROOT {
   static TClass *vectorlEvectorlEapvHitgRsPgR_Dictionary();
   static void vectorlEvectorlEapvHitgRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEapvHitgRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEapvHitgRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEapvHitgRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEapvHitgRsPgR(void *p);
   static void destruct_vectorlEvectorlEapvHitgRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<apvHit> >*)
   {
      vector<vector<apvHit> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<apvHit> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<apvHit> >", -2, "c++/v1/vector", 471,
                  typeid(vector<vector<apvHit> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEapvHitgRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<apvHit> >) );
      instance.SetNew(&new_vectorlEvectorlEapvHitgRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEapvHitgRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEapvHitgRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEapvHitgRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEapvHitgRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<apvHit> > >()));

      ::ROOT::AddClassAlternate("vector<vector<apvHit> >","std::__1::vector<std::__1::vector<apvHit, std::__1::allocator<apvHit> >, std::__1::allocator<std::__1::vector<apvHit, std::__1::allocator<apvHit> > > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vector<apvHit> >*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEapvHitgRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<apvHit> >*)nullptr)->GetClass();
      vectorlEvectorlEapvHitgRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEapvHitgRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEapvHitgRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<apvHit> > : new vector<vector<apvHit> >;
   }
   static void *newArray_vectorlEvectorlEapvHitgRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<apvHit> >[nElements] : new vector<vector<apvHit> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEapvHitgRsPgR(void *p) {
      delete ((vector<vector<apvHit> >*)p);
   }
   static void deleteArray_vectorlEvectorlEapvHitgRsPgR(void *p) {
      delete [] ((vector<vector<apvHit> >*)p);
   }
   static void destruct_vectorlEvectorlEapvHitgRsPgR(void *p) {
      typedef vector<vector<apvHit> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<apvHit> >

namespace ROOT {
   static TClass *vectorlEvectorlEapvClustergRsPgR_Dictionary();
   static void vectorlEvectorlEapvClustergRsPgR_TClassManip(TClass*);
   static void *new_vectorlEvectorlEapvClustergRsPgR(void *p = nullptr);
   static void *newArray_vectorlEvectorlEapvClustergRsPgR(Long_t size, void *p);
   static void delete_vectorlEvectorlEapvClustergRsPgR(void *p);
   static void deleteArray_vectorlEvectorlEapvClustergRsPgR(void *p);
   static void destruct_vectorlEvectorlEapvClustergRsPgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<vector<apvCluster> >*)
   {
      vector<vector<apvCluster> > *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<vector<apvCluster> >));
      static ::ROOT::TGenericClassInfo 
         instance("vector<vector<apvCluster> >", -2, "c++/v1/vector", 471,
                  typeid(vector<vector<apvCluster> >), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEvectorlEapvClustergRsPgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<vector<apvCluster> >) );
      instance.SetNew(&new_vectorlEvectorlEapvClustergRsPgR);
      instance.SetNewArray(&newArray_vectorlEvectorlEapvClustergRsPgR);
      instance.SetDelete(&delete_vectorlEvectorlEapvClustergRsPgR);
      instance.SetDeleteArray(&deleteArray_vectorlEvectorlEapvClustergRsPgR);
      instance.SetDestructor(&destruct_vectorlEvectorlEapvClustergRsPgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<vector<apvCluster> > >()));

      ::ROOT::AddClassAlternate("vector<vector<apvCluster> >","std::__1::vector<std::__1::vector<apvCluster, std::__1::allocator<apvCluster> >, std::__1::allocator<std::__1::vector<apvCluster, std::__1::allocator<apvCluster> > > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<vector<apvCluster> >*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEvectorlEapvClustergRsPgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<vector<apvCluster> >*)nullptr)->GetClass();
      vectorlEvectorlEapvClustergRsPgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEvectorlEapvClustergRsPgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEvectorlEapvClustergRsPgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<apvCluster> > : new vector<vector<apvCluster> >;
   }
   static void *newArray_vectorlEvectorlEapvClustergRsPgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<vector<apvCluster> >[nElements] : new vector<vector<apvCluster> >[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEvectorlEapvClustergRsPgR(void *p) {
      delete ((vector<vector<apvCluster> >*)p);
   }
   static void deleteArray_vectorlEvectorlEapvClustergRsPgR(void *p) {
      delete [] ((vector<vector<apvCluster> >*)p);
   }
   static void destruct_vectorlEvectorlEapvClustergRsPgR(void *p) {
      typedef vector<vector<apvCluster> > current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<vector<apvCluster> >

namespace ROOT {
   static TClass *vectorlEunsignedsPintgR_Dictionary();
   static void vectorlEunsignedsPintgR_TClassManip(TClass*);
   static void *new_vectorlEunsignedsPintgR(void *p = nullptr);
   static void *newArray_vectorlEunsignedsPintgR(Long_t size, void *p);
   static void delete_vectorlEunsignedsPintgR(void *p);
   static void deleteArray_vectorlEunsignedsPintgR(void *p);
   static void destruct_vectorlEunsignedsPintgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<unsigned int>*)
   {
      vector<unsigned int> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<unsigned int>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<unsigned int>", -2, "c++/v1/vector", 471,
                  typeid(vector<unsigned int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEunsignedsPintgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<unsigned int>) );
      instance.SetNew(&new_vectorlEunsignedsPintgR);
      instance.SetNewArray(&newArray_vectorlEunsignedsPintgR);
      instance.SetDelete(&delete_vectorlEunsignedsPintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEunsignedsPintgR);
      instance.SetDestructor(&destruct_vectorlEunsignedsPintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<unsigned int> >()));

      ::ROOT::AddClassAlternate("vector<unsigned int>","std::__1::vector<unsigned int, std::__1::allocator<unsigned int> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<unsigned int>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEunsignedsPintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<unsigned int>*)nullptr)->GetClass();
      vectorlEunsignedsPintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEunsignedsPintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEunsignedsPintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<unsigned int> : new vector<unsigned int>;
   }
   static void *newArray_vectorlEunsignedsPintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<unsigned int>[nElements] : new vector<unsigned int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEunsignedsPintgR(void *p) {
      delete ((vector<unsigned int>*)p);
   }
   static void deleteArray_vectorlEunsignedsPintgR(void *p) {
      delete [] ((vector<unsigned int>*)p);
   }
   static void destruct_vectorlEunsignedsPintgR(void *p) {
      typedef vector<unsigned int> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<unsigned int>

namespace ROOT {
   static TClass *vectorlEstringgR_Dictionary();
   static void vectorlEstringgR_TClassManip(TClass*);
   static void *new_vectorlEstringgR(void *p = nullptr);
   static void *newArray_vectorlEstringgR(Long_t size, void *p);
   static void delete_vectorlEstringgR(void *p);
   static void deleteArray_vectorlEstringgR(void *p);
   static void destruct_vectorlEstringgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<string>*)
   {
      vector<string> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<string>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<string>", -2, "c++/v1/vector", 471,
                  typeid(vector<string>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEstringgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<string>) );
      instance.SetNew(&new_vectorlEstringgR);
      instance.SetNewArray(&newArray_vectorlEstringgR);
      instance.SetDelete(&delete_vectorlEstringgR);
      instance.SetDeleteArray(&deleteArray_vectorlEstringgR);
      instance.SetDestructor(&destruct_vectorlEstringgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<string> >()));

      ::ROOT::AddClassAlternate("vector<string>","std::__1::vector<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char> >, std::__1::allocator<std::__1::basic_string<char, std::__1::char_traits<char>, std::__1::allocator<char> > > >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<string>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEstringgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<string>*)nullptr)->GetClass();
      vectorlEstringgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEstringgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEstringgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string> : new vector<string>;
   }
   static void *newArray_vectorlEstringgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<string>[nElements] : new vector<string>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEstringgR(void *p) {
      delete ((vector<string>*)p);
   }
   static void deleteArray_vectorlEstringgR(void *p) {
      delete [] ((vector<string>*)p);
   }
   static void destruct_vectorlEstringgR(void *p) {
      typedef vector<string> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<string>

namespace ROOT {
   static TClass *vectorlEshortgR_Dictionary();
   static void vectorlEshortgR_TClassManip(TClass*);
   static void *new_vectorlEshortgR(void *p = nullptr);
   static void *newArray_vectorlEshortgR(Long_t size, void *p);
   static void delete_vectorlEshortgR(void *p);
   static void deleteArray_vectorlEshortgR(void *p);
   static void destruct_vectorlEshortgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<short>*)
   {
      vector<short> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<short>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<short>", -2, "c++/v1/vector", 471,
                  typeid(vector<short>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEshortgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<short>) );
      instance.SetNew(&new_vectorlEshortgR);
      instance.SetNewArray(&newArray_vectorlEshortgR);
      instance.SetDelete(&delete_vectorlEshortgR);
      instance.SetDeleteArray(&deleteArray_vectorlEshortgR);
      instance.SetDestructor(&destruct_vectorlEshortgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<short> >()));

      ::ROOT::AddClassAlternate("vector<short>","std::__1::vector<short, std::__1::allocator<short> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<short>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEshortgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<short>*)nullptr)->GetClass();
      vectorlEshortgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEshortgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEshortgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<short> : new vector<short>;
   }
   static void *newArray_vectorlEshortgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<short>[nElements] : new vector<short>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEshortgR(void *p) {
      delete ((vector<short>*)p);
   }
   static void deleteArray_vectorlEshortgR(void *p) {
      delete [] ((vector<short>*)p);
   }
   static void destruct_vectorlEshortgR(void *p) {
      typedef vector<short> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<short>

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
         instance("vector<int>", -2, "c++/v1/vector", 471,
                  typeid(vector<int>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEintgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<int>) );
      instance.SetNew(&new_vectorlEintgR);
      instance.SetNewArray(&newArray_vectorlEintgR);
      instance.SetDelete(&delete_vectorlEintgR);
      instance.SetDeleteArray(&deleteArray_vectorlEintgR);
      instance.SetDestructor(&destruct_vectorlEintgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<int> >()));

      ::ROOT::AddClassAlternate("vector<int>","std::__1::vector<int, std::__1::allocator<int> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<int>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEintgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<int>*)nullptr)->GetClass();
      vectorlEintgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEintgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEintgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int> : new vector<int>;
   }
   static void *newArray_vectorlEintgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<int>[nElements] : new vector<int>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEintgR(void *p) {
      delete ((vector<int>*)p);
   }
   static void deleteArray_vectorlEintgR(void *p) {
      delete [] ((vector<int>*)p);
   }
   static void destruct_vectorlEintgR(void *p) {
      typedef vector<int> current_t;
      ((current_t*)p)->~current_t();
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
         instance("vector<float>", -2, "c++/v1/vector", 471,
                  typeid(vector<float>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEfloatgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<float>) );
      instance.SetNew(&new_vectorlEfloatgR);
      instance.SetNewArray(&newArray_vectorlEfloatgR);
      instance.SetDelete(&delete_vectorlEfloatgR);
      instance.SetDeleteArray(&deleteArray_vectorlEfloatgR);
      instance.SetDestructor(&destruct_vectorlEfloatgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<float> >()));

      ::ROOT::AddClassAlternate("vector<float>","std::__1::vector<float, std::__1::allocator<float> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<float>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEfloatgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<float>*)nullptr)->GetClass();
      vectorlEfloatgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEfloatgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEfloatgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float> : new vector<float>;
   }
   static void *newArray_vectorlEfloatgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<float>[nElements] : new vector<float>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEfloatgR(void *p) {
      delete ((vector<float>*)p);
   }
   static void deleteArray_vectorlEfloatgR(void *p) {
      delete [] ((vector<float>*)p);
   }
   static void destruct_vectorlEfloatgR(void *p) {
      typedef vector<float> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<float>

namespace ROOT {
   static TClass *vectorlEdoublegR_Dictionary();
   static void vectorlEdoublegR_TClassManip(TClass*);
   static void *new_vectorlEdoublegR(void *p = nullptr);
   static void *newArray_vectorlEdoublegR(Long_t size, void *p);
   static void delete_vectorlEdoublegR(void *p);
   static void deleteArray_vectorlEdoublegR(void *p);
   static void destruct_vectorlEdoublegR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<double>*)
   {
      vector<double> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<double>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<double>", -2, "c++/v1/vector", 471,
                  typeid(vector<double>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEdoublegR_Dictionary, isa_proxy, 4,
                  sizeof(vector<double>) );
      instance.SetNew(&new_vectorlEdoublegR);
      instance.SetNewArray(&newArray_vectorlEdoublegR);
      instance.SetDelete(&delete_vectorlEdoublegR);
      instance.SetDeleteArray(&deleteArray_vectorlEdoublegR);
      instance.SetDestructor(&destruct_vectorlEdoublegR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<double> >()));

      ::ROOT::AddClassAlternate("vector<double>","std::__1::vector<double, std::__1::allocator<double> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<double>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEdoublegR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<double>*)nullptr)->GetClass();
      vectorlEdoublegR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEdoublegR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEdoublegR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double> : new vector<double>;
   }
   static void *newArray_vectorlEdoublegR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<double>[nElements] : new vector<double>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEdoublegR(void *p) {
      delete ((vector<double>*)p);
   }
   static void deleteArray_vectorlEdoublegR(void *p) {
      delete [] ((vector<double>*)p);
   }
   static void destruct_vectorlEdoublegR(void *p) {
      typedef vector<double> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<double>

namespace ROOT {
   static TClass *vectorlEchargR_Dictionary();
   static void vectorlEchargR_TClassManip(TClass*);
   static void *new_vectorlEchargR(void *p = nullptr);
   static void *newArray_vectorlEchargR(Long_t size, void *p);
   static void delete_vectorlEchargR(void *p);
   static void deleteArray_vectorlEchargR(void *p);
   static void destruct_vectorlEchargR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<char>*)
   {
      vector<char> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<char>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<char>", -2, "c++/v1/vector", 471,
                  typeid(vector<char>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEchargR_Dictionary, isa_proxy, 4,
                  sizeof(vector<char>) );
      instance.SetNew(&new_vectorlEchargR);
      instance.SetNewArray(&newArray_vectorlEchargR);
      instance.SetDelete(&delete_vectorlEchargR);
      instance.SetDeleteArray(&deleteArray_vectorlEchargR);
      instance.SetDestructor(&destruct_vectorlEchargR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<char> >()));

      ::ROOT::AddClassAlternate("vector<char>","std::__1::vector<char, std::__1::allocator<char> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<char>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEchargR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<char>*)nullptr)->GetClass();
      vectorlEchargR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEchargR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEchargR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<char> : new vector<char>;
   }
   static void *newArray_vectorlEchargR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<char>[nElements] : new vector<char>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEchargR(void *p) {
      delete ((vector<char>*)p);
   }
   static void deleteArray_vectorlEchargR(void *p) {
      delete [] ((vector<char>*)p);
   }
   static void destruct_vectorlEchargR(void *p) {
      typedef vector<char> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<char>

namespace ROOT {
   static TClass *vectorlEapvTrackgR_Dictionary();
   static void vectorlEapvTrackgR_TClassManip(TClass*);
   static void *new_vectorlEapvTrackgR(void *p = nullptr);
   static void *newArray_vectorlEapvTrackgR(Long_t size, void *p);
   static void delete_vectorlEapvTrackgR(void *p);
   static void deleteArray_vectorlEapvTrackgR(void *p);
   static void destruct_vectorlEapvTrackgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<apvTrack>*)
   {
      vector<apvTrack> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<apvTrack>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<apvTrack>", -2, "c++/v1/vector", 471,
                  typeid(vector<apvTrack>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEapvTrackgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<apvTrack>) );
      instance.SetNew(&new_vectorlEapvTrackgR);
      instance.SetNewArray(&newArray_vectorlEapvTrackgR);
      instance.SetDelete(&delete_vectorlEapvTrackgR);
      instance.SetDeleteArray(&deleteArray_vectorlEapvTrackgR);
      instance.SetDestructor(&destruct_vectorlEapvTrackgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<apvTrack> >()));

      ::ROOT::AddClassAlternate("vector<apvTrack>","std::__1::vector<apvTrack, std::__1::allocator<apvTrack> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<apvTrack>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEapvTrackgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<apvTrack>*)nullptr)->GetClass();
      vectorlEapvTrackgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEapvTrackgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEapvTrackgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<apvTrack> : new vector<apvTrack>;
   }
   static void *newArray_vectorlEapvTrackgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<apvTrack>[nElements] : new vector<apvTrack>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEapvTrackgR(void *p) {
      delete ((vector<apvTrack>*)p);
   }
   static void deleteArray_vectorlEapvTrackgR(void *p) {
      delete [] ((vector<apvTrack>*)p);
   }
   static void destruct_vectorlEapvTrackgR(void *p) {
      typedef vector<apvTrack> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<apvTrack>

namespace ROOT {
   static TClass *vectorlEapvHitgR_Dictionary();
   static void vectorlEapvHitgR_TClassManip(TClass*);
   static void *new_vectorlEapvHitgR(void *p = nullptr);
   static void *newArray_vectorlEapvHitgR(Long_t size, void *p);
   static void delete_vectorlEapvHitgR(void *p);
   static void deleteArray_vectorlEapvHitgR(void *p);
   static void destruct_vectorlEapvHitgR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<apvHit>*)
   {
      vector<apvHit> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<apvHit>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<apvHit>", -2, "c++/v1/vector", 471,
                  typeid(vector<apvHit>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEapvHitgR_Dictionary, isa_proxy, 4,
                  sizeof(vector<apvHit>) );
      instance.SetNew(&new_vectorlEapvHitgR);
      instance.SetNewArray(&newArray_vectorlEapvHitgR);
      instance.SetDelete(&delete_vectorlEapvHitgR);
      instance.SetDeleteArray(&deleteArray_vectorlEapvHitgR);
      instance.SetDestructor(&destruct_vectorlEapvHitgR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<apvHit> >()));

      ::ROOT::AddClassAlternate("vector<apvHit>","std::__1::vector<apvHit, std::__1::allocator<apvHit> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<apvHit>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEapvHitgR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<apvHit>*)nullptr)->GetClass();
      vectorlEapvHitgR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEapvHitgR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEapvHitgR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<apvHit> : new vector<apvHit>;
   }
   static void *newArray_vectorlEapvHitgR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<apvHit>[nElements] : new vector<apvHit>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEapvHitgR(void *p) {
      delete ((vector<apvHit>*)p);
   }
   static void deleteArray_vectorlEapvHitgR(void *p) {
      delete [] ((vector<apvHit>*)p);
   }
   static void destruct_vectorlEapvHitgR(void *p) {
      typedef vector<apvHit> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<apvHit>

namespace ROOT {
   static TClass *vectorlEapvClustergR_Dictionary();
   static void vectorlEapvClustergR_TClassManip(TClass*);
   static void *new_vectorlEapvClustergR(void *p = nullptr);
   static void *newArray_vectorlEapvClustergR(Long_t size, void *p);
   static void delete_vectorlEapvClustergR(void *p);
   static void deleteArray_vectorlEapvClustergR(void *p);
   static void destruct_vectorlEapvClustergR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<apvCluster>*)
   {
      vector<apvCluster> *ptr = nullptr;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<apvCluster>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<apvCluster>", -2, "c++/v1/vector", 471,
                  typeid(vector<apvCluster>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEapvClustergR_Dictionary, isa_proxy, 4,
                  sizeof(vector<apvCluster>) );
      instance.SetNew(&new_vectorlEapvClustergR);
      instance.SetNewArray(&newArray_vectorlEapvClustergR);
      instance.SetDelete(&delete_vectorlEapvClustergR);
      instance.SetDeleteArray(&deleteArray_vectorlEapvClustergR);
      instance.SetDestructor(&destruct_vectorlEapvClustergR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<apvCluster> >()));

      ::ROOT::AddClassAlternate("vector<apvCluster>","std::__1::vector<apvCluster, std::__1::allocator<apvCluster> >");
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<apvCluster>*)nullptr); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEapvClustergR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<apvCluster>*)nullptr)->GetClass();
      vectorlEapvClustergR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEapvClustergR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEapvClustergR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<apvCluster> : new vector<apvCluster>;
   }
   static void *newArray_vectorlEapvClustergR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<apvCluster>[nElements] : new vector<apvCluster>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEapvClustergR(void *p) {
      delete ((vector<apvCluster>*)p);
   }
   static void deleteArray_vectorlEapvClustergR(void *p) {
      delete [] ((vector<apvCluster>*)p);
   }
   static void destruct_vectorlEapvClustergR(void *p) {
      typedef vector<apvCluster> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<apvCluster>

namespace {
  void TriggerDictionaryInitialization_link_Impl() {
    static const char* headers[] = {
"link.h",
nullptr
    };
    static const char* includePaths[] = {
"/Applications/root_v6.22.00/include/",
"/Volumes/Andrew's Seagate Drive/testBeamDATA/TB_CERN_JINR_2022/dsosnov_apv/TB_CERN_JINR_2022/code/",
nullptr
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "link dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_AutoLoading_Map;
namespace std{inline namespace __1{template <class _Tp> class __attribute__((annotate("$clingAutoload$iosfwd")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}}
struct __attribute__((annotate("$clingAutoload$apv_cluster.h")))  __attribute__((annotate("$clingAutoload$link.h")))  apvHit;
namespace std{inline namespace __1{template <class _CharT> struct __attribute__((annotate("$clingAutoload$__string")))  __attribute__((annotate("$clingAutoload$string")))  char_traits;
}}
class __attribute__((annotate("$clingAutoload$apv_cluster.h")))  __attribute__((annotate("$clingAutoload$link.h")))  apvCluster;
class __attribute__((annotate("$clingAutoload$apv_cluster.h")))  __attribute__((annotate("$clingAutoload$link.h")))  apvTrack;
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "link dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "link.h"

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("link",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_link_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_link_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_link() {
  TriggerDictionaryInitialization_link_Impl();
}
