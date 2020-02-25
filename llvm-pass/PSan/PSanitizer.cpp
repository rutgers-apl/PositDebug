/*
This pass instruments all posit computations to do shadow
execution in higer precision.  We need to handle two types of
variables - memory allocated and registers. Memory allocated variables
can be easily mapped to shadow memory using address as the key for
that variable.  Since registers can be passed to functions and
returned from functions, we need to simulate run time stack.  We need
a way to find where that register is allocated in stack. Unique
indexes for instructions help us to determine where that varaible is
allocated in shadow stack.
*/
				
#include "PSanitizer.h"
#include "llvm/IR/CallSite.h"  
#include "llvm/IR/ConstantFolder.h"
#include "llvm/ADT/SCCIterator.h" 
#include "llvm/ADT/StringExtras.h"
#include "llvm/ADT/StringRef.h"
#include <string>


//enum TYPE{MPFR, Posit64, Posit32, Posit16, Posit8};
enum TYPE{Double, Posit32, MPFR, Posit16, Posit8};

static cl::opt<int> Precision("fpsan-precision",
    cl::desc("default mpfr precision is initialized to 64"),
    cl::Hidden, cl::init(512));

static cl::opt<int> ENV("fpsan-with-type",
    cl::desc("shadow execution with mpfr"),
    cl::Hidden, cl::init(2));

void PSanitizer::addFunctionsToList(std::string FN) {
	std::ofstream myfile;
	myfile.open("functions.txt", std::ios::out|std::ios::app);
	if (myfile.is_open()){
		myfile <<FN;
		myfile << "\n";
		myfile.close();
	}
}

//check name of the function and check if it is in list of functions given by 
//developer and return true else false.
bool PSanitizer::isListedFunction(StringRef FN, std::string FileName) {
	std::ifstream infile(FileName);
	std::string line;
	while (std::getline(infile, line)) {
		if (FN.compare(line) == 0){
			return true;
		}
	}
	return false;
}

bool PSanitizer::isPositCall(Instruction *I){
  if(I == NULL)
    return false;
  if (CallInst *CI = dyn_cast<CallInst>(I)){
    Function *Callee = CI->getCalledFunction();
    if (Callee) {
      if(isListedFunction(Callee->getName(), "positFunc.txt")){
        return true;
      }
      if(Callee->getName().startswith("rapl_p32_force_load")){
        return true;
      }
    }
  }
  return false;
}

void PSanitizer::createGEP(Function *F, AllocaInst *Alloca){
  Function::iterator Fit = F->begin();
  BasicBlock &BB = *Fit;
  Instruction *I = dyn_cast<Instruction>(Alloca);
  Module *M = F->getParent();
  Instruction *Next = getNextInstruction(I, &BB);

  if(F->getName().startswith("__cons")) //avoid constructor
    return;
  IRBuilder<> IRB(Next);
  Instruction *End;
  for (auto &BB : *F) {
    for (auto &I : BB) {
      if (dyn_cast<ReturnInst>(&I)){
        End = &I;
      }
    }
  }
  IRBuilder<> IRBE(End);
  int index = 0;

  Type* VoidTy = Type::getVoidTy(M->getContext());
  for (auto &BB : *F) {
    for (auto &I : BB) {
      if (CallInst *CI = dyn_cast<CallInst>(&I)){
        Function *Callee = CI->getCalledFunction();
        if (Callee) {
          if(isListedFunction(Callee->getName(), "positFunc.txt")){
            Value *Indices[] = {ConstantInt::get(Type::getInt32Ty(M->getContext()), 0),
              ConstantInt::get(Type::getInt32Ty(M->getContext()), index)};
            Value *BOGEP = IRB.CreateGEP(Alloca, Indices);
            GEPMap.insert(std::pair<Instruction*, Value*>(&I, BOGEP));

            FuncInit = M->getOrInsertFunction("pd_init_mpfr", VoidTy, MPtrTy);
            IRB.CreateCall(FuncInit, {BOGEP});

            FuncInit = M->getOrInsertFunction("pd_clear_mpfr", VoidTy, MPtrTy);
            IRBE.CreateCall(FuncInit, {BOGEP});
            index++;
          }
        }
      }
      else if (AllocaInst *AI = dyn_cast<AllocaInst>(&I)){
        Value *Indices[] = {ConstantInt::get(Type::getInt32Ty(M->getContext()), 0),
          ConstantInt::get(Type::getInt32Ty(M->getContext()), index)};
        Value *BOGEP = IRB.CreateGEP(Alloca, Indices);
        GEPMap.insert(std::pair<Instruction*, Value*>(&I, BOGEP));

        FuncInit = M->getOrInsertFunction("pd_init_mpfr", VoidTy, MPtrTy);
        IRB.CreateCall(FuncInit, {BOGEP});

        FuncInit = M->getOrInsertFunction("pd_clear_mpfr", VoidTy, MPtrTy);
        IRBE.CreateCall(FuncInit, {BOGEP});
        index++;
      }
      else if(PHINode *PN = dyn_cast<PHINode>(&I)){
        bool flag = false;
        for (unsigned PI = 0, PE = PN->getNumIncomingValues(); PI != PE; ++PI) {
          Value *IncValue = PN->getIncomingValue(PI);

          if (IncValue == PN) continue; //TODO
          if(isPositCall(dyn_cast<Instruction>(IncValue))){
            flag = true;
          }
        }
        if(flag){
          Value *Indices[] = {ConstantInt::get(Type::getInt32Ty(M->getContext()), 0),
            ConstantInt::get(Type::getInt32Ty(M->getContext()), index)};
          Value *BOGEP = IRB.CreateGEP(Alloca, Indices);

          GEPMap.insert(std::pair<Instruction*, Value*>(&I, BOGEP));

          FuncInit = M->getOrInsertFunction("pd_init_mpfr", VoidTy, MPtrTy);
          IRB.CreateCall(FuncInit, {BOGEP});

          FuncInit = M->getOrInsertFunction("pd_clear_mpfr", VoidTy, MPtrTy);
          IRBE.CreateCall(FuncInit, {BOGEP});
          index++;
        }
      }
    }
      }
}

AllocaInst * PSanitizer::createAlloca(Function *F, size_t InsCount){
  Function::iterator Fit = F->begin();
  BasicBlock &BB = *Fit; 
  BasicBlock::iterator BBit = BB.begin();
  Instruction *First = &*BBit;
  IRBuilder<> IRB(First);
  Module *M = F->getParent();

  Instruction *End;
  for (auto &BB : *F) {	
    for (auto &I : BB) {
      if (dyn_cast<ReturnInst>(&I)){
        End = &I;
      }
    }
  }
  IRBuilder<> IRBE(End);

  //AllocaInst *Alloca = IRB.CreateAlloca(ArrayType::get(MPFRTy, InsCount),
  AllocaInst *Alloca = IRB.CreateAlloca(ArrayType::get(Real, InsCount),
      nullptr);
  return Alloca;
}

void PSanitizer::createMpfrAlloca(Function *F){
  long TotalArg = 1;
  long TotalAlloca = 0;
  for (Function::arg_iterator ait = F->arg_begin(), aend = F->arg_end();
      ait != aend; ++ait) {
    Argument *A = &*ait;
    ArgMap.insert(std::pair<Argument*, long>(A, TotalArg));
    TotalArg++;
  }

  FuncTotalArg.insert(std::pair<Function*, long>(F, TotalArg));
  TotalArg = 1;
  for (auto &BB : *F) {
    for (auto &I : BB) {
       if (CallInst *CI = dyn_cast<CallInst>(&I)){
        Function *Callee = CI->getCalledFunction();
        if (Callee) {
          if(Callee->getName().endswith("force_load")){
            TotalAlloca++;
          }
          else if(isListedFunction(Callee->getName(), "positFunc.txt")){
            TotalAlloca++;
          }
          else if(isListedFunction(Callee->getName(), "functions.txt")){
            TotalAlloca++;
          }
        }
      }
      else if(PHINode *PN = dyn_cast<PHINode>(&I)){
        bool flag = false;
        for (unsigned PI = 0, PE = PN->getNumIncomingValues(); PI != PE; ++PI) {
          Value *IncValue = PN->getIncomingValue(PI);
          if (IncValue == PN) continue; //TODO
          if(isPositCall(dyn_cast<Instruction>(IncValue))){
            TotalAlloca++;
            flag = true;
          }
        }
        if(flag) 
          TotalAlloca++;
      }
    }
  }
 // AllocaInst *Alloca = createAlloca(F, TotalAlloca);
//  createGEP(F, Alloca);
  TotalAlloca = 0;
}



//We need to call runtime function once floating point computation is executed, since
//we need to pass computed result to runtime to compare. To do that we need to instrument
//call after insted of before.
Instruction*
PSanitizer::getNextInstruction(Instruction *I, BasicBlock *BB){
  Instruction *Next;
  for (BasicBlock::iterator BBit = BB->begin(), BBend = BB->end();
       BBit != BBend; ++BBit) {
    Next = &*BBit;
    if(I == Next){
      Next = &*(++BBit);
      break;
    }
  }
  return Next;
}

Instruction*
PSanitizer::getNextInstructionNotPhi(Instruction *I, BasicBlock *BB){
  Instruction *Next;
    for (auto &I : *BB) {
      if(!isa<PHINode>(I)){
        Next = &I;
        break;
      }
  }
  return Next;
}

void PSanitizer::findInterestingFunctions(Function *F){

  bool flag = false;
  
  for (auto &BB : *F) {
    for (auto &I : BB) {
      if (CallInst *CI = dyn_cast<CallInst>(&I)){
	Function *Callee = CI->getCalledFunction();
	if (Callee) {
	  if(isListedFunction(Callee->getName(), "mathFunc.txt")){
	    flag = true;
	  }
	  if(isListedFunction(Callee->getName(), "positFunc.txt")){
	    flag = true;
	  }
	  if(Callee->getName().startswith("rapl_p32_set_arg")){
	    flag = true;
	  }
	}
      }
    }
  }

  if(isListedFunction(F->getName(), "wrapperFunc.txt")) return;
  if(flag){
    std::string name = F->getName();
    addFunctionsToList(name);
  }
}

void PSanitizer::handleFuncMainInit(Function *F){
  Function::iterator Fit = F->begin();
  BasicBlock &BB = *Fit;
  BasicBlock::iterator BBit = BB.begin();                                                                                                                                   
  Instruction *First = &*BBit;

  Module *M = F->getParent();
  IRBuilder<> IRB(First);

  Type* VoidTy = Type::getVoidTy(M->getContext());
  Type* Int64Ty = Type::getInt64Ty(M->getContext());

  Constant* Prec = ConstantInt::get(Type::getInt64Ty(M->getContext()), Precision);
  Finish = M->getOrInsertFunction("pd_init", VoidTy);
  long TotIns = 0;

  IRB.CreateCall(Finish, {});
}

void PSanitizer::handleMainRet(Instruction *I, Function *F){
  Module *M = F->getParent();
  IRBuilder<> IRB(I);
  Type* VoidTy = Type::getVoidTy(M->getContext());
  Finish = M->getOrInsertFunction("pd_finish", VoidTy);
  IRB.CreateCall(Finish, {});
}


void
PSanitizer::handleCallInst (CallInst *CI,
				BasicBlock *BB,
				Function *F,
				std::string CallName) {
  Instruction *I = dyn_cast<Instruction>(CI);
  Module *M = F->getParent();
  Instruction *Next = getNextInstruction(I, BB);
  IRBuilder<> IRBN(Next);
  IRBuilder<> IRB(I);

  Function *Callee = CI->getCalledFunction();

  Type* Int64Ty = Type::getInt64Ty(M->getContext());
  Type* VoidTy = Type::getVoidTy(M->getContext());
  Type* PtrVoidTy = PointerType::getUnqual(Type::getInt8Ty(M->getContext()));

  long InsIndex;
//  Value *BOGEP = GEPMap.at(CI);

  //if(isReturnPositType(CI)){
 //   FuncInit = M->getOrInsertFunction("__get_return", VoidTy, MPtrTy);
  //  IRBN.CreateCall(FuncInit, {BOGEP});
   // MInsMap.insert(std::pair<Instruction*, Instruction*>(I, dyn_cast<Instruction>(BOGEP)));
 // }

  size_t NumOperands = CI->getNumArgOperands();
  Value *Op[NumOperands];
  Type *OpTy[NumOperands];
  CallSite CS(I);
  for(int i = 0; i < NumOperands; i++){
    Op[i] = CI->getArgOperand(i);
    if (CS.paramHasAttr(i, Attribute::ByVal)){
      if(PointerType *Ty = dyn_cast<PointerType>(Op[i]->getType())){
        if (StructType *STy = dyn_cast<StructType>(Ty->getElementType())) {
          if(STy->getName() == "struct.posit_t"){
            Instruction *OpIns = dyn_cast<Instruction>(Op[i]);
            Value *OpIdx = handleOperand(I, Op[i], F);
            BitCastInst* BCToAddr = new BitCastInst(Op[i], 
                PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
            Constant* ArgNo = ConstantInt::get(Type::getInt64Ty(M->getContext()), i+1);
            AddFunArg = M->getOrInsertFunction("pd_set_argument", VoidTy, Int64Ty, PtrVoidTy);
            IRB.CreateCall(AddFunArg, {ArgNo, BCToAddr});
          }
        }
      }
    }
    if(!Op[i]->getType()->isPointerTy()){
      /*
        Instruction *OpIns = dyn_cast<Instruction>(Op[i]);
        Value *OpIdx = handleOperand(I, Op[i], F);
        Constant* ArgNo = ConstantInt::get(Type::getInt64Ty(M->getContext()), i+1);
        AddFunArg = M->getOrInsertFunction("__set_arg", VoidTy, Int64Ty, MPtrTy);
        IRB.CreateCall(AddFunArg, {ArgNo, OpIdx});
        */
    }
  }
}

void PSanitizer::handleFuncInit(Function *F){
  Function::iterator Fit = F->begin();
  BasicBlock &BB = *Fit; 
  BasicBlock::iterator BBit = BB.begin();
  Instruction *First = &*BBit;

  Module *M = F->getParent();
  IRBuilder<> IRB(First);
  Type* VoidTy = Type::getVoidTy(M->getContext());
  Type* Int64Ty = Type::getInt64Ty(M->getContext());

  FuncInit = M->getOrInsertFunction("pd_func_init", VoidTy, Int64Ty);
  long TotalArgs = FuncTotalArg.at(F);
  Constant* ConsTotIns = ConstantInt::get(Type::getInt64Ty(M->getContext()), TotalArgs); 
  IRB.CreateCall(FuncInit, {ConsTotIns});
}

ConstantInt*
PSanitizer::GetInstId(Instruction* I) {
  // Get unique instruction id
  MDNode* uniqueIdMDNode = I->getMetadata("psan_inst_id");
  if (uniqueIdMDNode == NULL) {
    llvm::errs()<<"ERROR: Cannot find instruction ID\n";
    exit(1);
  }
  Metadata* uniqueIdMetadata = uniqueIdMDNode->getOperand(0).get();
  ConstantAsMetadata* uniqueIdMD = dyn_cast<ConstantAsMetadata>(uniqueIdMetadata);
  Constant* uniqueIdConstant = uniqueIdMD->getValue();
  return dyn_cast<ConstantInt>(uniqueIdConstant);
}

void
PSanitizer::handlePositLibFunc (CallInst *CI,
				BasicBlock *BB,
				Function *F,
				Function * Callee) {
  Instruction *I = dyn_cast<Instruction>(CI);
  Instruction *Next = getNextInstruction(dyn_cast<Instruction>(CI), BB);
  IRBuilder<> IRB(Next);
  Module *M = F->getParent();
  
  Type* VoidTy = Type::getVoidTy(M->getContext());
  IntegerType* Int64Ty = Type::getInt64Ty(M->getContext());
  IntegerType* Int32Ty = Type::getInt32Ty(M->getContext());
  IntegerType* Int1Ty = Type::getInt1Ty(M->getContext());
  Type* PtrVoidTy = PointerType::getUnqual(Type::getInt8Ty(M->getContext()));

  // Get unique instruction id
  ConstantInt* instId = GetInstId(I);

  // Try to see if debug information is available
  const DebugLoc &instDebugLoc = I->getDebugLoc();
  bool debugInfoAvail = false;;
  unsigned int lineNum = 0;
  unsigned int colNum = 0;
  if (instDebugLoc) {
    debugInfoAvail = true;
    lineNum = instDebugLoc.getLine();
    colNum = instDebugLoc.getCol();
    if (lineNum == 0 && colNum == 0) debugInfoAvail = false;
  }

  ConstantInt* debugInfoAvailable = ConstantInt::get(Int1Ty, debugInfoAvail);
  ConstantInt* lineNumber = ConstantInt::get(Int32Ty, lineNum);
  ConstantInt* colNumber = ConstantInt::get(Int32Ty, colNum);


  // Now handle each posit functions
  if(Callee->getName().startswith("p32_expf")){
    BitCastInst* BCRes = new BitCastInst(CI->getOperand(0), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    BitCastInst* BCOp1 = new BitCastInst(CI->getOperand(1), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);

    FuncInit = M->getOrInsertFunction("pd_mpfr_exp", VoidTy, PtrVoidTy, PtrVoidTy, CI->getType(), Int64Ty, Int1Ty, Int32Ty, Int32Ty);
    IRB.CreateCall(FuncInit, {BCOp1, BCRes, CI, instId, debugInfoAvailable, lineNumber, colNumber});

    FuncInit = M->getOrInsertFunction("pd_check_error", Int32Ty, PtrVoidTy, CI->getOperand(0)->getType());
    IRB.CreateCall(FuncInit, {BCRes, CI->getOperand(0)});
  }
  else if(Callee->getName().startswith("p32_fabsf") || Callee->getName().startswith("p32_fabs")){
    BitCastInst* BCRes = new BitCastInst(CI->getOperand(0), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    BitCastInst* BCOp1 = new BitCastInst(CI->getOperand(1), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    FuncInit = M->getOrInsertFunction("pd_mpfr_fabs", VoidTy, PtrVoidTy, PtrVoidTy, 
        CI->getOperand(0)->getType(), Int64Ty, Int1Ty, Int32Ty, Int32Ty);
    IRB.CreateCall(FuncInit, {BCOp1, BCRes, CI->getOperand(0), instId, 
                  debugInfoAvailable, lineNumber, colNumber});

    FuncInit = M->getOrInsertFunction("pd_check_error", Int32Ty, PtrVoidTy, CI->getOperand(0)->getType());
    IRB.CreateCall(FuncInit, {BCRes, CI->getOperand(0)});
    return;
  }
  else if(Callee->getName().startswith("rapl_convertDoubleToP")){

    BitCastInst* BCToAddr = new BitCastInst(CI->getOperand(0), 
      PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    FuncInit = M->getOrInsertFunction("pd_set_const", VoidTy, PtrVoidTy, CI->getOperand(1)->getType(), Int32Ty);
    IRB.CreateCall(FuncInit, {BCToAddr, CI->getOperand(1), lineNumber});
    return;
  }
  else if(Callee->getName().endswith("_lt")){
    BitCastInst* BCOp1 = new BitCastInst(CI->getOperand(0), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    BitCastInst* BCOp2 = new BitCastInst(CI->getOperand(1), 
      PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    FuncInit = M->getOrInsertFunction("pd_mpfr_lt", Int1Ty, PtrVoidTy, PtrVoidTy, Int1Ty, Int64Ty, Int1Ty, Int32Ty, Int32Ty);
    Value *RealBr = IRB.CreateCall(FuncInit, {BCOp1, BCOp2, CI, instId, debugInfoAvailable, lineNumber, colNumber});

    FuncInit = M->getOrInsertFunction("pd_check_branch", Int32Ty, Int1Ty, Int1Ty, PtrVoidTy, PtrVoidTy);
    IRB.CreateCall(FuncInit, {RealBr, CI, BCOp1, BCOp2});
    return;
  }
  else if(Callee->getName().endswith("_le")){
    BitCastInst* BCOp1 = new BitCastInst(CI->getOperand(0), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    BitCastInst* BCOp2 = new BitCastInst(CI->getOperand(1), 
      PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    FuncInit = M->getOrInsertFunction("pd_mpfr_le", Int1Ty, PtrVoidTy, PtrVoidTy, Int1Ty, Int64Ty, Int1Ty, Int32Ty, Int32Ty);
    Value *RealBr = IRB.CreateCall(FuncInit, {BCOp1, BCOp2, CI, instId, debugInfoAvailable, lineNumber, colNumber});

    FuncInit = M->getOrInsertFunction("pd_check_branch", Int32Ty, Int1Ty, Int1Ty, PtrVoidTy, PtrVoidTy);
    IRB.CreateCall(FuncInit, {RealBr, CI, BCOp1, BCOp2});
    return;
  }
  else if(Callee->getName().endswith("_eq")){
    BitCastInst* BCOp1 = new BitCastInst(CI->getOperand(0), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    BitCastInst* BCOp2 = new BitCastInst(CI->getOperand(1), 
      PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    FuncInit = M->getOrInsertFunction("pd_mpfr_eq", Int1Ty, PtrVoidTy, PtrVoidTy, Int1Ty, Int64Ty, Int1Ty, Int32Ty, Int32Ty);
    Value *RealBr = IRB.CreateCall(FuncInit, {BCOp1, BCOp2, CI, instId, debugInfoAvailable, lineNumber, colNumber});

    FuncInit = M->getOrInsertFunction("pd_check_branch", Int32Ty, Int1Ty, Int1Ty, PtrVoidTy, PtrVoidTy);
    IRB.CreateCall(FuncInit, {RealBr, CI, BCOp1, BCOp2});
    return;
  }
  else if(Callee->getName().endswith("to_i32")){
    BitCastInst* BCOp1 = new BitCastInst(CI->getOperand(0), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    std::string Name = Callee->getName();
    FuncInit = M->getOrInsertFunction("pd_mpfr_"+Name, Int32Ty, PtrVoidTy, Int32Ty, Int64Ty, Int1Ty, Int32Ty, Int32Ty);
    Value *RealC = IRB.CreateCall(FuncInit, {BCOp1, CI, instId, debugInfoAvailable, lineNumber, colNumber});

    FuncInit = M->getOrInsertFunction("pd_check_conversion", Int32Ty, Int32Ty, Int32Ty, PtrVoidTy);
    IRB.CreateCall(FuncInit, {RealC, CI, BCOp1});
    return;
  }
  else if(Callee->getName().endswith("atan")){
    BitCastInst* BCRes = new BitCastInst(CI->getOperand(0), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    BitCastInst* BCOp1 = new BitCastInst(CI->getOperand(1), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    FuncInit = M->getOrInsertFunction("pd_mpfr_atan", VoidTy, PtrVoidTy, PtrVoidTy, 
        CI->getOperand(0)->getType(), Int64Ty, Int1Ty, Int32Ty, Int32Ty);
    IRB.CreateCall(FuncInit, {BCOp1, BCRes, CI->getOperand(0), instId, 
                  debugInfoAvailable, lineNumber, colNumber});

    FuncInit = M->getOrInsertFunction("pd_check_error", Int32Ty, PtrVoidTy, CI->getOperand(0)->getType());
    IRB.CreateCall(FuncInit, {BCRes, CI->getOperand(0)});
    return;
  }
  else if(Callee->getName().endswith("sin")){
    BitCastInst* BCRes = new BitCastInst(CI->getOperand(0), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    BitCastInst* BCOp1 = new BitCastInst(CI->getOperand(1), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    FuncInit = M->getOrInsertFunction("pd_mpfr_sin", VoidTy, PtrVoidTy, PtrVoidTy, 
        CI->getOperand(0)->getType(), Int64Ty, Int1Ty, Int32Ty, Int32Ty);
    IRB.CreateCall(FuncInit, {BCOp1, BCRes, CI->getOperand(0), instId, 
                  debugInfoAvailable, lineNumber, colNumber});

    FuncInit = M->getOrInsertFunction("pd_check_error", Int32Ty, PtrVoidTy, CI->getOperand(0)->getType());
    IRB.CreateCall(FuncInit, {BCRes, CI->getOperand(0)});
    return;
  }
  else if(Callee->getName().endswith("cos")){
    BitCastInst* BCRes = new BitCastInst(CI->getOperand(0), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    BitCastInst* BCOp1 = new BitCastInst(CI->getOperand(1), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    FuncInit = M->getOrInsertFunction("pd_mpfr_cos", VoidTy, PtrVoidTy, PtrVoidTy, 
        CI->getOperand(0)->getType(), Int64Ty, Int1Ty, Int32Ty, Int32Ty);
    IRB.CreateCall(FuncInit, {BCOp1, BCRes, CI->getOperand(0), instId, 
                  debugInfoAvailable, lineNumber, colNumber});

    FuncInit = M->getOrInsertFunction("pd_check_error", Int32Ty, PtrVoidTy, CI->getOperand(0)->getType());
    IRB.CreateCall(FuncInit, {BCRes, CI->getOperand(0)});
    return;
  }
  else if(Callee->getName().endswith("tan")){
    BitCastInst* BCRes = new BitCastInst(CI->getOperand(0), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    BitCastInst* BCOp1 = new BitCastInst(CI->getOperand(1), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    FuncInit = M->getOrInsertFunction("pd_mpfr_tan", VoidTy, PtrVoidTy, PtrVoidTy, 
        CI->getOperand(0)->getType(), Int64Ty, Int1Ty, Int32Ty, Int32Ty);
    IRB.CreateCall(FuncInit, {BCOp1, BCRes, CI->getOperand(0), instId, 
                  debugInfoAvailable, lineNumber, colNumber});

    FuncInit = M->getOrInsertFunction("pd_check_error", Int32Ty, PtrVoidTy, CI->getOperand(0)->getType());
    IRB.CreateCall(FuncInit, {BCRes, CI->getOperand(0)});
    return;
  }
  else if(Callee->getName().endswith("floor")){
    BitCastInst* BCRes = new BitCastInst(CI->getOperand(0), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    BitCastInst* BCOp1 = new BitCastInst(CI->getOperand(1), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    FuncInit = M->getOrInsertFunction("pd_mpfr_floor", VoidTy, PtrVoidTy, PtrVoidTy, 
        CI->getOperand(0)->getType(), Int64Ty, Int1Ty, Int32Ty, Int32Ty);
    IRB.CreateCall(FuncInit, {BCOp1, BCRes, CI->getOperand(0), instId, 
                  debugInfoAvailable, lineNumber, colNumber});

    FuncInit = M->getOrInsertFunction("pd_check_error", Int32Ty, PtrVoidTy, CI->getOperand(0)->getType());
    IRB.CreateCall(FuncInit, {BCRes, CI->getOperand(0)});
    return;
  }
  else if(Callee->getName().endswith("sqrt")){
    BitCastInst* BCRes = new BitCastInst(CI->getOperand(0), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    BitCastInst* BCOp1 = new BitCastInst(CI->getOperand(1), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    FuncInit = M->getOrInsertFunction("pd_mpfr_sqrt", VoidTy, PtrVoidTy, PtrVoidTy, 
        CI->getOperand(0)->getType(), Int64Ty, Int1Ty, Int32Ty, Int32Ty);
    IRB.CreateCall(FuncInit, {BCOp1, BCRes, CI->getOperand(0), instId, 
                  debugInfoAvailable, lineNumber, colNumber});

    FuncInit = M->getOrInsertFunction("pd_check_error", Int32Ty, PtrVoidTy, CI->getOperand(0)->getType());
    IRB.CreateCall(FuncInit, {BCRes, CI->getOperand(0)});
    return;
  }
  else if(Callee->getName().endswith("pow")){
    BitCastInst* BCRes = new BitCastInst(CI->getOperand(0), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    BitCastInst* BCOp1 = new BitCastInst(CI->getOperand(1), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    BitCastInst* BCOp2 = new BitCastInst(CI->getOperand(2), 
      PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    FuncInit = M->getOrInsertFunction("pd_mpfr_pow", VoidTy, PtrVoidTy, PtrVoidTy, PtrVoidTy, 
        CI->getType(), Int64Ty, Int1Ty, Int32Ty, Int32Ty);
    IRB.CreateCall(FuncInit, {BCOp1, BCOp2, BCRes, CI, instId, debugInfoAvailable, lineNumber, colNumber});

    FuncInit = M->getOrInsertFunction("pd_check_error", Int32Ty, PtrVoidTy, CI->getOperand(0)->getType());
    IRB.CreateCall(FuncInit, {BCRes, CI->getOperand(0)});
    return;
  }
  else{
    std::string Name = Callee->getName();

    Value *OP1 = CI->getOperand(1);
    Value *OP2 = CI->getOperand(2);

    BitCastInst* BCRes = new BitCastInst(CI->getOperand(0), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    BitCastInst* BCOp1 = new BitCastInst(CI->getOperand(1), 
        PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    BitCastInst* BCOp2 = new BitCastInst(CI->getOperand(2), 
      PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
    HandleFunc = M->getOrInsertFunction("pd_mpfr_"+Name, VoidTy, PtrVoidTy, PtrVoidTy, PtrVoidTy, 
        OP1->getType(), OP2->getType(), CI->getOperand(0)->getType(), Int64Ty, Int1Ty, Int32Ty, Int32Ty);

    IRB.CreateCall(HandleFunc, {BCOp1, BCOp2, BCRes, OP1, OP2, CI->getOperand(0), 
        instId, debugInfoAvailable, lineNumber, colNumber});

    FuncInit = M->getOrInsertFunction("pd_check_error", Int32Ty, PtrVoidTy, CI->getOperand(0)->getType());
    IRB.CreateCall(FuncInit, {BCRes, CI->getOperand(0)});
    return;
  }
}

void
PSanitizer::handleMathLibFunc (CallInst *CI,
				BasicBlock *BB,
				Function *F,
				std::string CallName) {
  Instruction *I = dyn_cast<Instruction>(CI);
  Instruction *Next = getNextInstruction(dyn_cast<Instruction>(CI), BB);
  IRBuilder<> IRB(Next);
  Module *M = F->getParent();
  
  Type* VoidTy = Type::getVoidTy(M->getContext());
  
  Type* Int64Ty = Type::getInt64Ty(M->getContext());
  Value *OP = CI->getOperand(0);
  
  Value *BOGEP = GEPMap.at(I);
  
  Type *OpTy = OP->getType();
  
  Value* ConsIdx1 = handleOperand(I, OP, F);

  HandleFunc = M->getOrInsertFunction("pd_mpfr_"+CallName, VoidTy, OpTy, MPtrTy, OpTy, MPtrTy);
  IRB.CreateCall(HandleFunc, {OP, ConsIdx1, CI, BOGEP});

  MInsMap.insert(std::pair<Value*, Instruction*>(CI->getOperand(0), dyn_cast<Instruction>(BOGEP)));
}

Value*
PSanitizer::handleOperand(Instruction *I, Value* OP, Function *F){

  Module *M = F->getParent();
  long Idx = 0;
	
  Value* ConsInsIndex;
  Instruction *OpIns = dyn_cast<Instruction>(OP);	
  //Either operand is 
  if(MInsMap.count(OP) != 0){
    ConsInsIndex = MInsMap.at(OP);
  }
  else if (GEPMap.count(OpIns) != 0){
    ConsInsIndex = GEPMap.at(OpIns);
  }
  else if(isa<Argument>(OP) && (MArgMap.count(dyn_cast<Argument>(OP)) != 0)){
    ConsInsIndex = MArgMap.at(dyn_cast<Argument>(OP));
  }
  else{
    //errs()<<"\nError !!! operand not found OP:"<<*OP<<"\n";
    ConsInsIndex = ConstantPointerNull::get(PointerType::get(Real, 0));
  }
  return ConsInsIndex;
}

void PSanitizer::handleRAPLSetReturn(CallInst *CI, BasicBlock *BB, Function *F){
  Instruction *I = dyn_cast<Instruction>(CI);
  Instruction *Next = getNextInstruction(I, BB);
  IRBuilder<> IRB(Next);
  Module *M = F->getParent();

  LLVMContext &C = F->getContext();
  Value *OP1 = CI->getOperand(0);

  Type* Int64Ty = Type::getInt64Ty(M->getContext());

  Type* VoidTy = Type::getVoidTy(M->getContext());
  Type* PtrVoidTy = PointerType::getUnqual(Type::getInt8Ty(M->getContext()));
  
  Instruction *OpIns = dyn_cast<Instruction>(OP1);
  //TODO: do we need to check for bitcast for store?
  Value *OP = CI->getOperand(0);
  Value *OpIdx = handleOperand(CI, OP1, F);
  long TotalArgs = FuncTotalArg.at(F);
  Constant* TotalArgsConst = ConstantInt::get(Type::getInt64Ty(M->getContext()), TotalArgs); 
  AddFunArg = M->getOrInsertFunction("pd_set_return", VoidTy, MPtrTy, Int64Ty, OP1->getType());
  IRB.CreateCall(AddFunArg, {OpIdx, TotalArgsConst, OP1});
}

void PSanitizer::handleRAPLGetReturn(CallInst *CI, BasicBlock *BB, Function *F){
  Instruction *I = dyn_cast<Instruction>(CI);
  Instruction *Next = getNextInstruction(I, BB);
  IRBuilder<> IRB(Next);
  Module *M = F->getParent();

  LLVMContext &C = F->getContext();
  Value *OP1 = CI->getOperand(0);

  Type* Int64Ty = Type::getInt64Ty(M->getContext());

  Type* VoidTy = Type::getVoidTy(M->getContext());
  Type* PtrVoidTy = PointerType::getUnqual(Type::getInt8Ty(M->getContext()));
  
  Instruction *OpIns = dyn_cast<Instruction>(OP1);
  //TODO: do we need to check for bitcast for store?
  Value *BOGEP = GEPMap.at(CI);
  FuncInit = M->getOrInsertFunction("pd_get_return", VoidTy, MPtrTy);
  IRB.CreateCall(FuncInit, {BOGEP});
  MInsMap.insert(std::pair<Value*, Instruction*>(CI->getOperand(0), dyn_cast<Instruction>(BOGEP)));
}

void PSanitizer::handleMemcpy(CallInst *CI, BasicBlock *BB, Function *F){
  Instruction *I = dyn_cast<Instruction>(CI);
  Instruction *Next = getNextInstruction(I, BB);
  IRBuilder<> IRB(Next);
  Module *M = F->getParent();

  Type* VoidTy = Type::getVoidTy(M->getContext());
  Type* PtrVoidTy = PointerType::getUnqual(Type::getInt8Ty(M->getContext()));
  BitCastInst* BCOp1 = new BitCastInst(CI->getOperand(1),
              PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);

  Value *OP1 = CI->getOperand(0);
  Value *OP2 = CI->getOperand(1);
//  OP1->stripPointerCasts()->dump();
 // OP1->getType()->dump();
  //if (BitCastInst *BI = dyn_cast<BitCastInst>(OP1)){
//  StructType *STy = cast<StructType>(OP1->stripPointerCasts()->getType()->getPointerElementType());
  if (BitCastInst *BI = dyn_cast<BitCastInst>(CI->getOperand(0))){
    Type *BITy = BI->getOperand(0)->getType();
    if(BITy->getPointerElementType()->getTypeID() == Type::StructTyID){
      StructType *STy = cast<StructType>(BITy->getPointerElementType());
      if(STy->getName() == "struct.posit_t"){
        AddFunArg = M->getOrInsertFunction("pd_handle_memcpy", VoidTy, PtrVoidTy, PtrVoidTy);
        IRB.CreateCall(AddFunArg, {OP1, OP2});
      }
    }
  }
  else if (BitCastInst *BI = dyn_cast<BitCastInst>(CI->getOperand(1))){
    Type *BITy = BI->getOperand(0)->getType();
    if(BITy->getPointerElementType()->getTypeID() == Type::StructTyID){
      StructType *STy = cast<StructType>(BITy->getPointerElementType());
      if(STy->getName() == "struct.posit_t"){
        AddFunArg = M->getOrInsertFunction("pd_handle_memcpy", VoidTy, PtrVoidTy, PtrVoidTy);
        IRB.CreateCall(AddFunArg, {OP1, OP2});
      }
    }
  }
}

void PSanitizer::handleRAPLSetArg(CallInst *CI, BasicBlock *BB, Function *F){
  Instruction *I = dyn_cast<Instruction>(CI);
  Instruction *Next = getNextInstruction(I, BB);
  IRBuilder<> IRB(Next);
  Module *M = F->getParent();

  LLVMContext &C = F->getContext();
  Value *OP1 = CI->getOperand(0);
  Value *OP2 = CI->getOperand(1);

  Type* Int32Ty = Type::getInt32Ty(M->getContext());

  Type* VoidTy = Type::getVoidTy(M->getContext());
  Type* PtrVoidTy = PointerType::getUnqual(Type::getInt8Ty(M->getContext()));
  

  Instruction *OpIns = dyn_cast<Instruction>(OP1);
  //TODO: do we need to check for bitcast for store?

  Value *OpIdx = handleOperand(I, OP1, F);
  AddFunArg = M->getOrInsertFunction("pd_set_argument", VoidTy, Int32Ty, MPtrTy);
  IRB.CreateCall(AddFunArg, {OP2, OpIdx});
}

void PSanitizer::handleRAPLGetArg(Function *F){
  Function::iterator Fit = F->begin();
  BasicBlock &BB = *Fit; 
  BasicBlock::iterator BBit = BB.begin();
  Instruction *First = &*BBit;
  IRBuilder<> IRB(First);
  Module *M = F->getParent();

  Type* VoidTy = Type::getVoidTy(M->getContext());
  Type* Int64Ty = Type::getInt64Ty(M->getContext());
  Type* PtrVoidTy = PointerType::getUnqual(Type::getInt8Ty(M->getContext()));

  for (Function::arg_iterator ait = F->arg_begin(), aend = F->arg_end();
      ait != aend; ++ait) {
    Argument *A = &*ait;
    if (A->hasByValAttr()){
      if(PointerType *Ty = dyn_cast<PointerType>(A->getType())){
        if (StructType *STy = dyn_cast<StructType>(Ty->getElementType())) {
          if(STy->getName() == "struct.posit_t"){
            //TODO: do we need to check for bitcast for store?
            size_t Idx =  ArgMap.at(A);

            long TotalArgs = FuncTotalArg.at(F);
            BitCastInst* BCToAddr = new BitCastInst(A,
                PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", First);

            Constant* ArgNo = ConstantInt::get(Type::getInt64Ty(M->getContext()), TotalArgs-Idx);
            SetRealTemp = M->getOrInsertFunction("pd_get_argument", VoidTy, Int64Ty, PtrVoidTy);
            Value* ConsInsIndex = IRB.CreateCall(SetRealTemp, {ArgNo, BCToAddr});
            MArgMap.insert(std::pair<Argument*, Instruction*>(A, dyn_cast<Instruction>(ConsInsIndex)));
          }
        }
      }
    }
  }
}

void PSanitizer::handlePhi(PHINode *PN, BasicBlock *BB, Function *F){
  Module *M = F->getParent();
  Type* Int64Ty = Type::getInt64Ty(M->getContext());
  Type* VoidTy = Type::getVoidTy(M->getContext());
  IRBuilder<> IRB(dyn_cast<Instruction>(dyn_cast<Instruction>(PN)));

  PHINode* iPHI = IRB.CreatePHI (MPtrTy, 2);
  MInsMap.insert(std::pair<Value*, Instruction*>(PN, iPHI));
  NewPhiMap.insert(std::pair<Instruction*, Instruction*>(dyn_cast<Instruction>(PN), iPHI));

  Instruction* Next;
  Next = getNextInstructionNotPhi(PN, BB);
  IRBuilder<> IRBN(Next);
  if(GEPMap.count(PN) != 0){
    Value *BOGEP = GEPMap.at(PN);

    AddFunArg = M->getOrInsertFunction("pd_copy_phi", VoidTy, MPtrTy, MPtrTy);
    IRBN.CreateCall(AddFunArg, {iPHI, BOGEP});
  }
}

void PSanitizer::handleRAPLStore(CallInst *CI, BasicBlock *BB, Function *F){
  Instruction *I = dyn_cast<Instruction>(CI);
  Instruction *Next = getNextInstruction(I, BB);
  IRBuilder<> IRB(Next);
  Module *M = F->getParent();

  LLVMContext &C = F->getContext();
  Value *Addr = CI->getOperand(0);
  Value *OP = CI->getOperand(1);

  Type* Int64Ty = Type::getInt64Ty(M->getContext());

  Type* VoidTy = Type::getVoidTy(M->getContext());
  Type* PtrVoidTy = PointerType::getUnqual(Type::getInt8Ty(M->getContext()));
  
  Instruction *OpIns = dyn_cast<Instruction>(OP);
  //TODO: do we need to check for bitcast for store?
  BitCastInst* BCToAddr = new BitCastInst(Addr, 
      PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);
  Value* InsIndex = handleOperand(I, OP, F);
  SetRealTemp = M->getOrInsertFunction("pd_store_real", VoidTy, PtrVoidTy, MPtrTy);
  IRB.CreateCall(SetRealTemp, {BCToAddr, InsIndex});
}


void PSanitizer::handleNewPhi(Function *F){
  Module *M = F->getParent();
  Instruction* Next;
  long NumPhi = 0;
  BasicBlock *IBB, *BB;
  for(auto it = NewPhiMap.begin(); it != NewPhiMap.end(); ++it)
  {
    if(PHINode *PN = dyn_cast<PHINode>(it->first)){
      PHINode* iPHI = dyn_cast<PHINode>(it->second);
      for (unsigned PI = 0, PE = PN->getNumIncomingValues(); PI != PE; ++PI) {
        IBB = PN->getIncomingBlock(PI);
        Value *IncValue = PN->getIncomingValue(PI);
        BB = PN->getParent();

        if (IncValue == PN) continue; //TODO
        Value* InsIndex = handleOperand(it->first, IncValue, F);
        iPHI->addIncoming(InsIndex, IBB);
      }
    }
  }
}

//handleReturn should call set_return before mpfr_clear
void PSanitizer::handleReturn(ReturnInst *RI, BasicBlock *BB, Function *F){
  Instruction *Ins = dyn_cast<Instruction>(RI);
  Module *M = F->getParent();

  Type* Int64Ty = Type::getInt64Ty(M->getContext());
  Type* VoidTy = Type::getVoidTy(M->getContext());

  Value* OpIdx;
  //Find first mpfr clear
  for (auto &BB : *F) {
    for (auto &I : BB) {
      if (CallInst *CI = dyn_cast<CallInst>(&I)){
        Function *Callee = CI->getCalledFunction();
        if(Callee && Callee->getName() == "pd_clear_mpfr"){
            Ins = &I;
          break;
        }
      }
    }
  }
  IRBuilder<> IRB(Ins);
  FuncInit = M->getOrInsertFunction("pd_func_exit", VoidTy, Int64Ty);
  long TotalArgs = FuncTotalArg.at(F);
  Constant* ConsTotIns = ConstantInt::get(Type::getInt64Ty(M->getContext()), TotalArgs); 
  IRB.CreateCall(FuncInit, {ConsTotIns});
}

void PSanitizer::handleRAPLLoad(CallInst *CI, BasicBlock *BB, Function *F){
  Instruction *I = dyn_cast<Instruction>(CI);
  Module *M = F->getParent();
  Instruction *Next = getNextInstruction(I, BB);
  IRBuilder<> IRB(Next);

  LLVMContext &C = F->getContext();

  Type* PtrVoidTy = PointerType::getUnqual(Type::getInt8Ty(C));
  Type* VoidTy = Type::getVoidTy(C);
  Type* Int1Ty = Type::getInt1Ty(C);
  Type* Int64Ty = Type::getInt64Ty(C);

  Value *Addr = CI->getOperand(0);
  Value *BOGEP = GEPMap.at(CI);

  BitCastInst* BCToAddr = new BitCastInst(Addr, 
      PointerType::getUnqual(Type::getInt8Ty(M->getContext())),"", I);

  LoadCall = M->getOrInsertFunction("pd_load_shadow", VoidTy, MPtrTy, PtrVoidTy);
  Value* LoadI = IRB.CreateCall(LoadCall, {BOGEP, BCToAddr});
  MInsMap.insert(std::pair<Value*, Instruction*>(CI->getOperand(0), dyn_cast<Instruction>(BOGEP)));
}

void PSanitizer::handleIns(Instruction *I, BasicBlock *BB, Function *F){
  //instrument interesting instructions

  if (CallInst *CI = dyn_cast<CallInst>(I)){
    Function *Callee = CI->getCalledFunction();
    if (!Callee) return;

    if(Callee->getName().startswith("llvm.memcpy")){
      handleMemcpy(CI, BB, F);
    }
    else if(Callee->getName().endswith("set_arg")){
      handleRAPLSetArg(CI, BB, F);
    }
    else if(Callee->getName().endswith("get_ret")){
      handleRAPLGetReturn(CI, BB, F);
    }
    else if(Callee->getName().endswith("set_ret")){
      handleRAPLSetReturn(CI, BB, F);
    }
    else if(isListedFunction(Callee->getName(), "positFunc.txt")){
      handlePositLibFunc(CI, BB, F, Callee);
    }
    else if(isListedFunction(Callee->getName(), "mathFunc.txt")) {
      handleMathLibFunc(CI, BB, F, Callee->getName());
    }
    else if(isListedFunction(Callee->getName(), "functions.txt")) {
      if(!Callee->getName().startswith("rapl_")){
         handleCallInst(CI, BB, F, Callee->getName());
      }
    }
  }
}

bool PSanitizer::runOnModule(Module &M) {
  LLVMContext &C = M.getContext();

  // Create MPFR type and MPFR pointer type
  if(ENV == MPFR){
    StructType* MPFRTy1 = StructType::create(M.getContext(), "struct.fpsan_mpfr");
    MPFRTy1->setBody({Type::getInt64Ty(M.getContext()), Type::getInt32Ty(M.getContext()), 
                    Type::getInt64Ty(M.getContext()), Type::getInt64PtrTy(M.getContext())});

    MPFRTy = StructType::create(M.getContext(), "struct.f_mpfr");
    MPFRTy->setBody(llvm::ArrayType::get(MPFRTy1, 1));
//    MPtrTy = MPFRTy->getPointerTo();
    
    Real = StructType::create(M.getContext(), "struct.fpsan_real");
    RealPtr = Real->getPointerTo();
    Real->setBody({MPFRTy, Type::getDoubleTy(M.getContext()), Type::getInt64Ty(M.getContext()),
                  Type::getInt64Ty(M.getContext()), Type::getInt64Ty(M.getContext()), 
                  Type::getInt64Ty(M.getContext()), Type::getInt64Ty(M.getContext()),
                  Type::getInt64Ty(M.getContext()), Type::getInt64Ty(M.getContext()),
                  Type::getInt32Ty(M.getContext()),
                  RealPtr, RealPtr, Type::getInt64Ty(M.getContext()), 
                  Type::getInt64Ty(M.getContext()), Type::getInt1Ty(M.getContext())});

    MPtrTy = Real->getPointerTo();
  }
  
  //TODO::Iterate over global arrays to initialize shadow memory
  for (Module::global_iterator GVI = M.global_begin(), E = M.global_end();
               GVI != E; ) {
    GlobalVariable *GV = &*GVI++;
    if(GV->hasInitializer()){
      Constant *Init = GV->getInitializer();
    }
  }
  
  // Find functions that perform posit computation. No instrumentation
  // if the function does not perform any posit computations.
  for (auto &F : M) {
    if (F.isDeclaration()) continue;
    findInterestingFunctions(&F);
  }

  for (auto &F : M) {
    if (F.isDeclaration()) continue;
    if (!isListedFunction(F.getName(), "functions.txt")) continue;
    //All instrumented functions are listed in AllFuncList	
    AllFuncList.push_back(&F);
  }

  //instrument interesting instructions
  Instruction *LastPhi = NULL;
  int instId = 0;
  for (Function *F : reverse(AllFuncList)) {
    // Create MPFR metadata on the stack
    createMpfrAlloca(F);
    handleRAPLGetArg(F);

    //if argument is used in any floating point computation, then we
    //need to retrieve that argument from shadow stack.  Instead of
    //call __get_arg everytime opearnd is used, it is better to call
    //once in start of the function and remember the address of shadow
    //stack.

    if(F->getName() != "main"){
      //add func_init and func_exit in the start and end of the
      //function to set shadow stack variables
      if(!F->getName().startswith("__cons")){
        if(!F->getName().startswith("pd_init")){ //avoid constructor
          if(!F->getName().startswith("p16_")){
            if(!F->getName().startswith("p32_")){
              handleFuncInit(F);
            }
          }
        }
      }
    }
    
    for (auto &BB : *F) {
      for (auto &I : BB) {
        // For each instruction, add a unique id to the metadata of I
        LLVMContext& instContext = I.getContext();
        ConstantInt* instUniqueId = ConstantInt::get(Type::getInt64Ty(M.getContext()), instId);
        ConstantAsMetadata* uniqueId = ConstantAsMetadata::get(instUniqueId);
        MDNode* md = MDNode::get(instContext, uniqueId);
        I.setMetadata("psan_inst_id", md);
        instId++;

        if(PHINode *PN = dyn_cast<PHINode>(&I)){
          bool flag = false;
          for (unsigned PI = 0, PE = PN->getNumIncomingValues(); PI != PE; ++PI) {
            Value *IncValue = PN->getIncomingValue(PI);

            if (IncValue == PN) continue; //TODO
            if(isPositCall(dyn_cast<Instruction>(IncValue))){
              flag = true;
            }
          }
          if(flag){
            handlePhi(PN, &BB, F);
            LastPhi = &I;
          }
        }
        handleIns(&I, &BB, F);
      }
    }
    
    for (auto &BB : *F) {
      for (auto &I : BB) {
        if (ReturnInst *RI = dyn_cast<ReturnInst>(&I)){
          if(F->getName() != "main"){
            if(!F->getName().startswith("__cons"))
              if(!F->getName().startswith("__init")) //avoid constructor
                if(!F->getName().startswith("p32"))
                  handleReturn(RI, &BB, F);
          }
        }
      }
    }
    handleNewPhi(F);
    NewPhiMap.clear(); 
    MInsMap.clear(); 
    GEPMap.clear(); 
    ConsMap.clear(); 
  }
  
  for (auto &F : M) {
    if (F.isDeclaration()) continue;
    if(F.getName() == "main"){
      //add init and finish func in the start and end of the main
      //function to initialize shadow memory
      handleFuncMainInit(&F);
    }
    for (auto &BB : F) {
      for (auto &I : BB) {
        if (dyn_cast<ReturnInst>(&I)){
          if(F.getName() == "main"){
            handleMainRet(&I, &F);
          }
        }
      }
    }
  } 
 
  return true;
}



void addFPPass(const PassManagerBuilder &Builder, legacy::PassManagerBase &PM) {
  PM.add(new PSanitizer());
}

RegisterStandardPasses SOpt(PassManagerBuilder::EP_OptimizerLast,
			    addFPPass);
RegisterStandardPasses S(PassManagerBuilder::EP_EnabledOnOptLevel0,
                         addFPPass);

char PSanitizer::ID = 0;
static const RegisterPass<PSanitizer> Y("psan", "instrument fp operations", false, false);
