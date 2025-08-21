# mlwe
Rust implementation of a simplified version of module learning with errors (M-LWE), a lattice based post quantum cryptosystem.  It can be used for encryption, key-exchange, or digital signatures.

The repository currently implements encryption and decryption of byte strings.  There is a library that exports encrypt/decrypt functionality, with a test case to show that decryption recovers encrypted data.  There is also a binary that  that takes strings on the command line, which are encrypted then decrypted.
