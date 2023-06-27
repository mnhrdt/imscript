int cpu_new(float*, int,int,int);           // open cpu with float image data
void cpu_update(int, float*, int,int,int);  // update cpu window data
void cpu_send_key(int, int);                // send a key to cpu (by utf8 code)
void cpu_close(int);                        // close cpu window
