// Do NOT change. Changes will be lost next time file is generated

#define R__DICTIONARY_FILENAME MuEtreeDict
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

// The generated code does not explicitly qualifies STL entities
namespace std {} using namespace std;

// Header files passed as explicit arguments
#include "MuEtree.h"
#include "MuEana.h"
#include "MuEtree.h"
#include "MuEana.h"

// Header files passed via #pragma extra_include

namespace ROOT {
   static void *new_MuEcLcLMCpara(void *p = 0);
   static void *newArray_MuEcLcLMCpara(Long_t size, void *p);
   static void delete_MuEcLcLMCpara(void *p);
   static void deleteArray_MuEcLcLMCpara(void *p);
   static void destruct_MuEcLcLMCpara(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MuE::MCpara*)
   {
      ::MuE::MCpara *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MuE::MCpara >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MuE::MCpara", ::MuE::MCpara::Class_Version(), "MuEtree.h", 19,
                  typeid(::MuE::MCpara), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MuE::MCpara::Dictionary, isa_proxy, 4,
                  sizeof(::MuE::MCpara) );
      instance.SetNew(&new_MuEcLcLMCpara);
      instance.SetNewArray(&newArray_MuEcLcLMCpara);
      instance.SetDelete(&delete_MuEcLcLMCpara);
      instance.SetDeleteArray(&deleteArray_MuEcLcLMCpara);
      instance.SetDestructor(&destruct_MuEcLcLMCpara);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MuE::MCpara*)
   {
      return GenerateInitInstanceLocal((::MuE::MCpara*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::MuE::MCpara*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_MuEcLcLMCstat(void *p = 0);
   static void *newArray_MuEcLcLMCstat(Long_t size, void *p);
   static void delete_MuEcLcLMCstat(void *p);
   static void deleteArray_MuEcLcLMCstat(void *p);
   static void destruct_MuEcLcLMCstat(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MuE::MCstat*)
   {
      ::MuE::MCstat *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MuE::MCstat >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MuE::MCstat", ::MuE::MCstat::Class_Version(), "MuEtree.h", 51,
                  typeid(::MuE::MCstat), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MuE::MCstat::Dictionary, isa_proxy, 4,
                  sizeof(::MuE::MCstat) );
      instance.SetNew(&new_MuEcLcLMCstat);
      instance.SetNewArray(&newArray_MuEcLcLMCstat);
      instance.SetDelete(&delete_MuEcLcLMCstat);
      instance.SetDeleteArray(&deleteArray_MuEcLcLMCstat);
      instance.SetDestructor(&destruct_MuEcLcLMCstat);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MuE::MCstat*)
   {
      return GenerateInitInstanceLocal((::MuE::MCstat*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::MuE::MCstat*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_MuEcLcLSetup(void *p = 0);
   static void *newArray_MuEcLcLSetup(Long_t size, void *p);
   static void delete_MuEcLcLSetup(void *p);
   static void deleteArray_MuEcLcLSetup(void *p);
   static void destruct_MuEcLcLSetup(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MuE::Setup*)
   {
      ::MuE::Setup *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MuE::Setup >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MuE::Setup", ::MuE::Setup::Class_Version(), "MuEtree.h", 82,
                  typeid(::MuE::Setup), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MuE::Setup::Dictionary, isa_proxy, 4,
                  sizeof(::MuE::Setup) );
      instance.SetNew(&new_MuEcLcLSetup);
      instance.SetNewArray(&newArray_MuEcLcLSetup);
      instance.SetDelete(&delete_MuEcLcLSetup);
      instance.SetDeleteArray(&deleteArray_MuEcLcLSetup);
      instance.SetDestructor(&destruct_MuEcLcLSetup);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MuE::Setup*)
   {
      return GenerateInitInstanceLocal((::MuE::Setup*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::MuE::Setup*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_MuEcLcLP4(void *p = 0);
   static void *newArray_MuEcLcLP4(Long_t size, void *p);
   static void delete_MuEcLcLP4(void *p);
   static void deleteArray_MuEcLcLP4(void *p);
   static void destruct_MuEcLcLP4(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MuE::P4*)
   {
      ::MuE::P4 *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MuE::P4 >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MuE::P4", ::MuE::P4::Class_Version(), "MuEtree.h", 118,
                  typeid(::MuE::P4), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MuE::P4::Dictionary, isa_proxy, 4,
                  sizeof(::MuE::P4) );
      instance.SetNew(&new_MuEcLcLP4);
      instance.SetNewArray(&newArray_MuEcLcLP4);
      instance.SetDelete(&delete_MuEcLcLP4);
      instance.SetDeleteArray(&deleteArray_MuEcLcLP4);
      instance.SetDestructor(&destruct_MuEcLcLP4);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MuE::P4*)
   {
      return GenerateInitInstanceLocal((::MuE::P4*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::MuE::P4*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_MuEcLcLEvent(void *p = 0);
   static void *newArray_MuEcLcLEvent(Long_t size, void *p);
   static void delete_MuEcLcLEvent(void *p);
   static void deleteArray_MuEcLcLEvent(void *p);
   static void destruct_MuEcLcLEvent(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MuE::Event*)
   {
      ::MuE::Event *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MuE::Event >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MuE::Event", ::MuE::Event::Class_Version(), "MuEtree.h", 135,
                  typeid(::MuE::Event), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MuE::Event::Dictionary, isa_proxy, 4,
                  sizeof(::MuE::Event) );
      instance.SetNew(&new_MuEcLcLEvent);
      instance.SetNewArray(&newArray_MuEcLcLEvent);
      instance.SetDelete(&delete_MuEcLcLEvent);
      instance.SetDeleteArray(&deleteArray_MuEcLcLEvent);
      instance.SetDestructor(&destruct_MuEcLcLEvent);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MuE::Event*)
   {
      return GenerateInitInstanceLocal((::MuE::Event*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::MuE::Event*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_MuEcLcLKineVars(void *p = 0);
   static void *newArray_MuEcLcLKineVars(Long_t size, void *p);
   static void delete_MuEcLcLKineVars(void *p);
   static void deleteArray_MuEcLcLKineVars(void *p);
   static void destruct_MuEcLcLKineVars(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MuE::KineVars*)
   {
      ::MuE::KineVars *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MuE::KineVars >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MuE::KineVars", ::MuE::KineVars::Class_Version(), "MuEana.h", 14,
                  typeid(::MuE::KineVars), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MuE::KineVars::Dictionary, isa_proxy, 4,
                  sizeof(::MuE::KineVars) );
      instance.SetNew(&new_MuEcLcLKineVars);
      instance.SetNewArray(&newArray_MuEcLcLKineVars);
      instance.SetDelete(&delete_MuEcLcLKineVars);
      instance.SetDeleteArray(&deleteArray_MuEcLcLKineVars);
      instance.SetDestructor(&destruct_MuEcLcLKineVars);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MuE::KineVars*)
   {
      return GenerateInitInstanceLocal((::MuE::KineVars*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::MuE::KineVars*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_MuEcLcLPhoton(void *p = 0);
   static void *newArray_MuEcLcLPhoton(Long_t size, void *p);
   static void delete_MuEcLcLPhoton(void *p);
   static void deleteArray_MuEcLcLPhoton(void *p);
   static void destruct_MuEcLcLPhoton(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MuE::Photon*)
   {
      ::MuE::Photon *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MuE::Photon >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MuE::Photon", ::MuE::Photon::Class_Version(), "MuEana.h", 101,
                  typeid(::MuE::Photon), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MuE::Photon::Dictionary, isa_proxy, 4,
                  sizeof(::MuE::Photon) );
      instance.SetNew(&new_MuEcLcLPhoton);
      instance.SetNewArray(&newArray_MuEcLcLPhoton);
      instance.SetDelete(&delete_MuEcLcLPhoton);
      instance.SetDeleteArray(&deleteArray_MuEcLcLPhoton);
      instance.SetDestructor(&destruct_MuEcLcLPhoton);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MuE::Photon*)
   {
      return GenerateInitInstanceLocal((::MuE::Photon*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::MuE::Photon*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace ROOT {
   static void *new_MuEcLcLMuEana(void *p = 0);
   static void *newArray_MuEcLcLMuEana(Long_t size, void *p);
   static void delete_MuEcLcLMuEana(void *p);
   static void deleteArray_MuEcLcLMuEana(void *p);
   static void destruct_MuEcLcLMuEana(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const ::MuE::MuEana*)
   {
      ::MuE::MuEana *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TInstrumentedIsAProxy< ::MuE::MuEana >(0);
      static ::ROOT::TGenericClassInfo 
         instance("MuE::MuEana", ::MuE::MuEana::Class_Version(), "MuEana.h", 122,
                  typeid(::MuE::MuEana), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &::MuE::MuEana::Dictionary, isa_proxy, 4,
                  sizeof(::MuE::MuEana) );
      instance.SetNew(&new_MuEcLcLMuEana);
      instance.SetNewArray(&newArray_MuEcLcLMuEana);
      instance.SetDelete(&delete_MuEcLcLMuEana);
      instance.SetDeleteArray(&deleteArray_MuEcLcLMuEana);
      instance.SetDestructor(&destruct_MuEcLcLMuEana);
      return &instance;
   }
   TGenericClassInfo *GenerateInitInstance(const ::MuE::MuEana*)
   {
      return GenerateInitInstanceLocal((::MuE::MuEana*)0);
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const ::MuE::MuEana*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));
} // end of namespace ROOT

namespace MuE {
//______________________________________________________________________________
atomic_TClass_ptr MCpara::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MCpara::Class_Name()
{
   return "MuE::MCpara";
}

//______________________________________________________________________________
const char *MCpara::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MuE::MCpara*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MCpara::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MuE::MCpara*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MCpara::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MuE::MCpara*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MCpara::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MuE::MCpara*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace MuE
namespace MuE {
//______________________________________________________________________________
atomic_TClass_ptr MCstat::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MCstat::Class_Name()
{
   return "MuE::MCstat";
}

//______________________________________________________________________________
const char *MCstat::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MuE::MCstat*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MCstat::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MuE::MCstat*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MCstat::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MuE::MCstat*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MCstat::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MuE::MCstat*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace MuE
namespace MuE {
//______________________________________________________________________________
atomic_TClass_ptr Setup::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Setup::Class_Name()
{
   return "MuE::Setup";
}

//______________________________________________________________________________
const char *Setup::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MuE::Setup*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Setup::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MuE::Setup*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Setup::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MuE::Setup*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Setup::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MuE::Setup*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace MuE
namespace MuE {
//______________________________________________________________________________
atomic_TClass_ptr P4::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *P4::Class_Name()
{
   return "MuE::P4";
}

//______________________________________________________________________________
const char *P4::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MuE::P4*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int P4::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MuE::P4*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *P4::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MuE::P4*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *P4::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MuE::P4*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace MuE
namespace MuE {
//______________________________________________________________________________
atomic_TClass_ptr Event::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Event::Class_Name()
{
   return "MuE::Event";
}

//______________________________________________________________________________
const char *Event::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MuE::Event*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Event::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MuE::Event*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Event::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MuE::Event*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Event::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MuE::Event*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace MuE
namespace MuE {
//______________________________________________________________________________
atomic_TClass_ptr KineVars::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *KineVars::Class_Name()
{
   return "MuE::KineVars";
}

//______________________________________________________________________________
const char *KineVars::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MuE::KineVars*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int KineVars::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MuE::KineVars*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *KineVars::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MuE::KineVars*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *KineVars::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MuE::KineVars*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace MuE
namespace MuE {
//______________________________________________________________________________
atomic_TClass_ptr Photon::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *Photon::Class_Name()
{
   return "MuE::Photon";
}

//______________________________________________________________________________
const char *Photon::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MuE::Photon*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int Photon::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MuE::Photon*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *Photon::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MuE::Photon*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *Photon::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MuE::Photon*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace MuE
namespace MuE {
//______________________________________________________________________________
atomic_TClass_ptr MuEana::fgIsA(0);  // static to hold class pointer

//______________________________________________________________________________
const char *MuEana::Class_Name()
{
   return "MuE::MuEana";
}

//______________________________________________________________________________
const char *MuEana::ImplFileName()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MuE::MuEana*)0x0)->GetImplFileName();
}

//______________________________________________________________________________
int MuEana::ImplFileLine()
{
   return ::ROOT::GenerateInitInstanceLocal((const ::MuE::MuEana*)0x0)->GetImplFileLine();
}

//______________________________________________________________________________
TClass *MuEana::Dictionary()
{
   fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MuE::MuEana*)0x0)->GetClass();
   return fgIsA;
}

//______________________________________________________________________________
TClass *MuEana::Class()
{
   if (!fgIsA.load()) { R__LOCKGUARD(gInterpreterMutex); fgIsA = ::ROOT::GenerateInitInstanceLocal((const ::MuE::MuEana*)0x0)->GetClass(); }
   return fgIsA;
}

} // namespace MuE
namespace MuE {
//______________________________________________________________________________
void MCpara::Streamer(TBuffer &R__b)
{
   // Stream an object of class MuE::MCpara.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(MuE::MCpara::Class(),this);
   } else {
      R__b.WriteClassBuffer(MuE::MCpara::Class(),this);
   }
}

} // namespace MuE
namespace ROOT {
   // Wrappers around operator new
   static void *new_MuEcLcLMCpara(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::MuE::MCpara : new ::MuE::MCpara;
   }
   static void *newArray_MuEcLcLMCpara(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::MuE::MCpara[nElements] : new ::MuE::MCpara[nElements];
   }
   // Wrapper around operator delete
   static void delete_MuEcLcLMCpara(void *p) {
      delete ((::MuE::MCpara*)p);
   }
   static void deleteArray_MuEcLcLMCpara(void *p) {
      delete [] ((::MuE::MCpara*)p);
   }
   static void destruct_MuEcLcLMCpara(void *p) {
      typedef ::MuE::MCpara current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MuE::MCpara

namespace MuE {
//______________________________________________________________________________
void MCstat::Streamer(TBuffer &R__b)
{
   // Stream an object of class MuE::MCstat.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(MuE::MCstat::Class(),this);
   } else {
      R__b.WriteClassBuffer(MuE::MCstat::Class(),this);
   }
}

} // namespace MuE
namespace ROOT {
   // Wrappers around operator new
   static void *new_MuEcLcLMCstat(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::MuE::MCstat : new ::MuE::MCstat;
   }
   static void *newArray_MuEcLcLMCstat(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::MuE::MCstat[nElements] : new ::MuE::MCstat[nElements];
   }
   // Wrapper around operator delete
   static void delete_MuEcLcLMCstat(void *p) {
      delete ((::MuE::MCstat*)p);
   }
   static void deleteArray_MuEcLcLMCstat(void *p) {
      delete [] ((::MuE::MCstat*)p);
   }
   static void destruct_MuEcLcLMCstat(void *p) {
      typedef ::MuE::MCstat current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MuE::MCstat

namespace MuE {
//______________________________________________________________________________
void Setup::Streamer(TBuffer &R__b)
{
   // Stream an object of class MuE::Setup.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(MuE::Setup::Class(),this);
   } else {
      R__b.WriteClassBuffer(MuE::Setup::Class(),this);
   }
}

} // namespace MuE
namespace ROOT {
   // Wrappers around operator new
   static void *new_MuEcLcLSetup(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::MuE::Setup : new ::MuE::Setup;
   }
   static void *newArray_MuEcLcLSetup(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::MuE::Setup[nElements] : new ::MuE::Setup[nElements];
   }
   // Wrapper around operator delete
   static void delete_MuEcLcLSetup(void *p) {
      delete ((::MuE::Setup*)p);
   }
   static void deleteArray_MuEcLcLSetup(void *p) {
      delete [] ((::MuE::Setup*)p);
   }
   static void destruct_MuEcLcLSetup(void *p) {
      typedef ::MuE::Setup current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MuE::Setup

namespace MuE {
//______________________________________________________________________________
void P4::Streamer(TBuffer &R__b)
{
   // Stream an object of class MuE::P4.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(MuE::P4::Class(),this);
   } else {
      R__b.WriteClassBuffer(MuE::P4::Class(),this);
   }
}

} // namespace MuE
namespace ROOT {
   // Wrappers around operator new
   static void *new_MuEcLcLP4(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::MuE::P4 : new ::MuE::P4;
   }
   static void *newArray_MuEcLcLP4(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::MuE::P4[nElements] : new ::MuE::P4[nElements];
   }
   // Wrapper around operator delete
   static void delete_MuEcLcLP4(void *p) {
      delete ((::MuE::P4*)p);
   }
   static void deleteArray_MuEcLcLP4(void *p) {
      delete [] ((::MuE::P4*)p);
   }
   static void destruct_MuEcLcLP4(void *p) {
      typedef ::MuE::P4 current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MuE::P4

namespace MuE {
//______________________________________________________________________________
void Event::Streamer(TBuffer &R__b)
{
   // Stream an object of class MuE::Event.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(MuE::Event::Class(),this);
   } else {
      R__b.WriteClassBuffer(MuE::Event::Class(),this);
   }
}

} // namespace MuE
namespace ROOT {
   // Wrappers around operator new
   static void *new_MuEcLcLEvent(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::MuE::Event : new ::MuE::Event;
   }
   static void *newArray_MuEcLcLEvent(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::MuE::Event[nElements] : new ::MuE::Event[nElements];
   }
   // Wrapper around operator delete
   static void delete_MuEcLcLEvent(void *p) {
      delete ((::MuE::Event*)p);
   }
   static void deleteArray_MuEcLcLEvent(void *p) {
      delete [] ((::MuE::Event*)p);
   }
   static void destruct_MuEcLcLEvent(void *p) {
      typedef ::MuE::Event current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MuE::Event

namespace MuE {
//______________________________________________________________________________
void KineVars::Streamer(TBuffer &R__b)
{
   // Stream an object of class MuE::KineVars.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(MuE::KineVars::Class(),this);
   } else {
      R__b.WriteClassBuffer(MuE::KineVars::Class(),this);
   }
}

} // namespace MuE
namespace ROOT {
   // Wrappers around operator new
   static void *new_MuEcLcLKineVars(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::MuE::KineVars : new ::MuE::KineVars;
   }
   static void *newArray_MuEcLcLKineVars(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::MuE::KineVars[nElements] : new ::MuE::KineVars[nElements];
   }
   // Wrapper around operator delete
   static void delete_MuEcLcLKineVars(void *p) {
      delete ((::MuE::KineVars*)p);
   }
   static void deleteArray_MuEcLcLKineVars(void *p) {
      delete [] ((::MuE::KineVars*)p);
   }
   static void destruct_MuEcLcLKineVars(void *p) {
      typedef ::MuE::KineVars current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MuE::KineVars

namespace MuE {
//______________________________________________________________________________
void Photon::Streamer(TBuffer &R__b)
{
   // Stream an object of class MuE::Photon.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(MuE::Photon::Class(),this);
   } else {
      R__b.WriteClassBuffer(MuE::Photon::Class(),this);
   }
}

} // namespace MuE
namespace ROOT {
   // Wrappers around operator new
   static void *new_MuEcLcLPhoton(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::MuE::Photon : new ::MuE::Photon;
   }
   static void *newArray_MuEcLcLPhoton(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::MuE::Photon[nElements] : new ::MuE::Photon[nElements];
   }
   // Wrapper around operator delete
   static void delete_MuEcLcLPhoton(void *p) {
      delete ((::MuE::Photon*)p);
   }
   static void deleteArray_MuEcLcLPhoton(void *p) {
      delete [] ((::MuE::Photon*)p);
   }
   static void destruct_MuEcLcLPhoton(void *p) {
      typedef ::MuE::Photon current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MuE::Photon

namespace MuE {
//______________________________________________________________________________
void MuEana::Streamer(TBuffer &R__b)
{
   // Stream an object of class MuE::MuEana.

   if (R__b.IsReading()) {
      R__b.ReadClassBuffer(MuE::MuEana::Class(),this);
   } else {
      R__b.WriteClassBuffer(MuE::MuEana::Class(),this);
   }
}

} // namespace MuE
namespace ROOT {
   // Wrappers around operator new
   static void *new_MuEcLcLMuEana(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::MuE::MuEana : new ::MuE::MuEana;
   }
   static void *newArray_MuEcLcLMuEana(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) ::MuE::MuEana[nElements] : new ::MuE::MuEana[nElements];
   }
   // Wrapper around operator delete
   static void delete_MuEcLcLMuEana(void *p) {
      delete ((::MuE::MuEana*)p);
   }
   static void deleteArray_MuEcLcLMuEana(void *p) {
      delete [] ((::MuE::MuEana*)p);
   }
   static void destruct_MuEcLcLMuEana(void *p) {
      typedef ::MuE::MuEana current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class ::MuE::MuEana

namespace ROOT {
   static TClass *vectorlEMuEcLcLP4gR_Dictionary();
   static void vectorlEMuEcLcLP4gR_TClassManip(TClass*);
   static void *new_vectorlEMuEcLcLP4gR(void *p = 0);
   static void *newArray_vectorlEMuEcLcLP4gR(Long_t size, void *p);
   static void delete_vectorlEMuEcLcLP4gR(void *p);
   static void deleteArray_vectorlEMuEcLcLP4gR(void *p);
   static void destruct_vectorlEMuEcLcLP4gR(void *p);

   // Function generating the singleton type initializer
   static TGenericClassInfo *GenerateInitInstanceLocal(const vector<MuE::P4>*)
   {
      vector<MuE::P4> *ptr = 0;
      static ::TVirtualIsAProxy* isa_proxy = new ::TIsAProxy(typeid(vector<MuE::P4>));
      static ::ROOT::TGenericClassInfo 
         instance("vector<MuE::P4>", -2, "vector", 339,
                  typeid(vector<MuE::P4>), ::ROOT::Internal::DefineBehavior(ptr, ptr),
                  &vectorlEMuEcLcLP4gR_Dictionary, isa_proxy, 4,
                  sizeof(vector<MuE::P4>) );
      instance.SetNew(&new_vectorlEMuEcLcLP4gR);
      instance.SetNewArray(&newArray_vectorlEMuEcLcLP4gR);
      instance.SetDelete(&delete_vectorlEMuEcLcLP4gR);
      instance.SetDeleteArray(&deleteArray_vectorlEMuEcLcLP4gR);
      instance.SetDestructor(&destruct_vectorlEMuEcLcLP4gR);
      instance.AdoptCollectionProxyInfo(TCollectionProxyInfo::Generate(TCollectionProxyInfo::Pushback< vector<MuE::P4> >()));
      return &instance;
   }
   // Static variable to force the class initialization
   static ::ROOT::TGenericClassInfo *_R__UNIQUE_DICT_(Init) = GenerateInitInstanceLocal((const vector<MuE::P4>*)0x0); R__UseDummy(_R__UNIQUE_DICT_(Init));

   // Dictionary for non-ClassDef classes
   static TClass *vectorlEMuEcLcLP4gR_Dictionary() {
      TClass* theClass =::ROOT::GenerateInitInstanceLocal((const vector<MuE::P4>*)0x0)->GetClass();
      vectorlEMuEcLcLP4gR_TClassManip(theClass);
   return theClass;
   }

   static void vectorlEMuEcLcLP4gR_TClassManip(TClass* ){
   }

} // end of namespace ROOT

namespace ROOT {
   // Wrappers around operator new
   static void *new_vectorlEMuEcLcLP4gR(void *p) {
      return  p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<MuE::P4> : new vector<MuE::P4>;
   }
   static void *newArray_vectorlEMuEcLcLP4gR(Long_t nElements, void *p) {
      return p ? ::new((::ROOT::Internal::TOperatorNewHelper*)p) vector<MuE::P4>[nElements] : new vector<MuE::P4>[nElements];
   }
   // Wrapper around operator delete
   static void delete_vectorlEMuEcLcLP4gR(void *p) {
      delete ((vector<MuE::P4>*)p);
   }
   static void deleteArray_vectorlEMuEcLcLP4gR(void *p) {
      delete [] ((vector<MuE::P4>*)p);
   }
   static void destruct_vectorlEMuEcLcLP4gR(void *p) {
      typedef vector<MuE::P4> current_t;
      ((current_t*)p)->~current_t();
   }
} // end of namespace ROOT for class vector<MuE::P4>

namespace {
  void TriggerDictionaryInitialization_MuEtreeDict_Impl() {
    static const char* headers[] = {
"MuEtree.h",
"MuEana.h",
0
    };
    static const char* includePaths[] = {
"../code",
"/cvmfs/sft.cern.ch/lcg/releases/ROOT/v6.20.02-10e75/x86_64-centos7-gcc8-opt/include/",
"/home/LHCB-T3/espedicato/RobeTesi/provaECAL/test/",
0
    };
    static const char* fwdDeclCode = R"DICTFWDDCLS(
#line 1 "MuEtreeDict dictionary forward declarations' payload"
#pragma clang diagnostic ignored "-Wkeyword-compat"
#pragma clang diagnostic ignored "-Wignored-attributes"
#pragma clang diagnostic ignored "-Wreturn-type-c-linkage"
extern int __Cling_Autoloading_Map;
namespace MuE{class __attribute__((annotate("$clingAutoload$MuEtree.h")))  P4;}
namespace std{template <typename _Tp> class __attribute__((annotate("$clingAutoload$bits/allocator.h")))  __attribute__((annotate("$clingAutoload$string")))  allocator;
}
namespace MuE{class __attribute__((annotate("$clingAutoload$MuEtree.h")))  MCpara;}
namespace MuE{class __attribute__((annotate("$clingAutoload$MuEtree.h")))  MCstat;}
namespace MuE{class __attribute__((annotate("$clingAutoload$MuEtree.h")))  Setup;}
namespace MuE{class __attribute__((annotate("$clingAutoload$MuEtree.h")))  Event;}
namespace MuE{class __attribute__((annotate("$clingAutoload$MuEana.h")))  KineVars;}
namespace MuE{class __attribute__((annotate("$clingAutoload$MuEana.h")))  Photon;}
namespace MuE{class __attribute__((annotate("$clingAutoload$MuEana.h")))  MuEana;}
)DICTFWDDCLS";
    static const char* payloadCode = R"DICTPAYLOAD(
#line 1 "MuEtreeDict dictionary payload"


#define _BACKWARD_BACKWARD_WARNING_H
// Inline headers
#include "MuEtree.h"
#include "MuEana.h"
#include "MuEtree.h"
#include "MuEana.h"

#if defined(__CLING__) || defined(__CINT__)
#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class MuE::P4+;
#pragma link C++ class std::vector<MuE::P4>+;
#pragma link C++ class MuE::Event+;
#pragma link C++ class MuE::MCpara+;
#pragma link C++ class MuE::MCstat+;
#pragma link C++ class MuE::Setup+;
#pragma link C++ class MuE::KineVars+;
#pragma link C++ class MuE::Photon+;
#pragma link C++ class MuE::MuEana+;

#endif

#undef  _BACKWARD_BACKWARD_WARNING_H
)DICTPAYLOAD";
    static const char* classesHeaders[] = {
"MuE::Event", payloadCode, "@",
"MuE::KineVars", payloadCode, "@",
"MuE::MCpara", payloadCode, "@",
"MuE::MCstat", payloadCode, "@",
"MuE::MuEana", payloadCode, "@",
"MuE::P4", payloadCode, "@",
"MuE::Photon", payloadCode, "@",
"MuE::Setup", payloadCode, "@",
nullptr
};
    static bool isInitialized = false;
    if (!isInitialized) {
      TROOT::RegisterModule("MuEtreeDict",
        headers, includePaths, payloadCode, fwdDeclCode,
        TriggerDictionaryInitialization_MuEtreeDict_Impl, {}, classesHeaders, /*hasCxxModule*/false);
      isInitialized = true;
    }
  }
  static struct DictInit {
    DictInit() {
      TriggerDictionaryInitialization_MuEtreeDict_Impl();
    }
  } __TheDictionaryInitializer;
}
void TriggerDictionaryInitialization_MuEtreeDict() {
  TriggerDictionaryInitialization_MuEtreeDict_Impl();
}
