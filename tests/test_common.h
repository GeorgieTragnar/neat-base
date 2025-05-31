#pragma once

// Enable access to protected and private members for testing
// This must be defined BEFORE any project headers are included
#ifdef TESTING
#define protected public
#define private public
#endif

// Note: Include this header AFTER all standard library headers 
// but BEFORE any project headers in your test files