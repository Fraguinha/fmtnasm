; example hello world program

                global    _main                                                 ; entry point
                extern    _puts                                                 ; link with puts

                section   .data                                                 ; initialized data, constants

message:        db        "Hello, world!", 0                                    ; string to print, 10 -> '\n', 0 -> '\0'

                section   .text                                                 ; code

_main:
                push      rbp                                                   ; save rbp
                mov       rbp, rsp                                              ; create frame

                lea       rdi, [rel message]                                    ; 1st argument
                call      _puts                                                 ; puts(message)

                leave                                                           ; destroy frame
                ret                                                             ; return
