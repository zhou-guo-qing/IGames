﻿#ifndef iGameMacro_h
#define iGameMacro_h

#include <sstream>
#include "iGameFileUtils.h"

 // 命名空间的定义
#define IGAME_NAMESPACE iGame
#define IGAME_NAMESPACE_BEGIN namespace iGame {
#define IGAME_NAMESPACE_END }

#define I_OBJECT(TypeName)                         \
  TypeName(const TypeName &) = delete;             \
  TypeName & operator=(const TypeName &) = delete; \
  TypeName(TypeName &&) = delete;                  \
  TypeName & operator=(TypeName &&) = delete;      \
  using Self = TypeName;                           \
  using Pointer = SmartPointer<Self>;              \
  using ConstPointer = SmartPointer<const Self>; 

#define I_TEMPLATED_OBJECT(TypeName,TempName)      \
  TypeName(const TypeName &) = delete;             \
  TypeName & operator=(const TypeName &) = delete; \
  TypeName(TypeName &&) = delete;                  \
  TypeName & operator=(TypeName &&) = delete;      \
  using Self = TypeName<TempName>;                 \
  using Pointer = SmartPointer<Self>;              \
  using ConstPointer = SmartPointer<const Self>; 


#define igDebug(...) iGameDebugMacro(this, __VA_ARGS__)
#define igError(...) iGameErrorMacro(this, __VA_ARGS__)

#define iGameDebugMacro(self, msg)                                \
do {                                                              \
    std::stringstream _ss;                                        \
    _ss << msg;                                                   \
    std::string const _filename = GetFilenameName(__FILE__);      \
    OutputWindowDisplayDebugText(_filename.c_str(), __LINE__, _ss.str().c_str()); \
} while(0)

#define iGameErrorMacro(self, msg)                                \
do {                                                              \
    std::stringstream _ss;                                        \
    _ss << msg;                                                   \
    std::string const _filename = GetFilenameName(__FILE__);      \
    OutputWindowDisplayErrorText(_filename.c_str(), __LINE__, _ss.str().c_str()); \
    exit(0);                                                      \
} while(0)





#endif