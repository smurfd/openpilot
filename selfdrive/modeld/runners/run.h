#pragma once

#include "selfdrive/modeld/runners/runmodel.h"
#include "selfdrive/modeld/runners/snpemodel.h"

#if defined(QCOM) || defined(QCOM2)
#include "selfdrive/modeld/runners/thneedmodel.h"
#define DefaultRunModel SNPEModel
#else
#ifdef USE_ONNX_MODEL
#include "selfdrive/modeld/runners/onnxmodel.h"
#define DefaultRunModel ONNXModel
#else
#define DefaultRunModel SNPEModel
#endif
#endif
