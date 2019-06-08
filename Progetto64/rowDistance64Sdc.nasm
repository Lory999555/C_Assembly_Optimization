%include "sseutils64.nasm"

section .data
	
subb equ 16

dim equ 4


section .bss
    
section .text

global rowDistance64SdcA

rowDistance64SdcA:

    push		rbp				; salva il Base Pointer
    mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
    pushaq						; salva i registri generali

    mov r11,[rbp+subb]
    ;vprintreg r11
    
    ;mov dword [distance],0              ;distance=0
      
    imul rcx,dim                        ;prendo j
    imul rcx,r11                        ;j*sub
    imul rcx,r9                        ;j*k*subb

    imul rdx,dim                        ;prendo i
    imul rdx,r11                        ;i*subb
    add rdx,rcx                         ;j*k*subb+i*subb
    add rdx,rdi                         ;c[j*k*subb+i*subb]
 
    imul r8,dim                         ;prendo j_d
    imul r8,r11                        ;j_d*subb
    add r8,rcx                         ;j*k*subb+j_d*subb
    add r8,rdi                         ;c[j*k*subb+j_d*subb]
    
    ;printreg rdx
    ;printreg r9
    vmovaps ymm1,[rdx]                   ;c[j*k*sub+i*sub+0]
    vsubps ymm1,[r8]                    ;c[j*k*sub+i*sub+0] - c[j*k*sub+j_d*sub+0]
    vmulps ymm1,ymm1                     ;(..)^2
    ;printregyps ymm1
    ;mov rdi,dim
    sub r11,16                          ;subb-16

    mov r12,8                           ;z=8
    vxorps ymm2, ymm2
ciclo:
    cmp r12,r11
    jg resto

    vmovaps ymm0,[rdx+4*r12]             ;c[j*k*sub+i*sub+z]
    ;printregyps ymm0
    vmovaps ymm7,[rdx+4*r12+32]             ;c[j*k*sub+i*sub+z]
    ;printregyps ymm0
    vsubps ymm0,[r8+4*r12]              ;c[j*k*sub+i*sub+z] - c[j*k*sub+j_d*sub+z]
    vsubps ymm7,[r8+4*r12+32]              ;c[j*k*sub+i*sub+z] - c[j*k*sub+j_d*sub+z]
    ;printregyps ymm0
    vmulps ymm0,ymm0                     ;(..)^2
    vmulps ymm7,ymm7                     ;(..)^2
    ;printregyps ymm0
    vaddps ymm1,ymm0                     ;distance+=(..)^2
    vaddps ymm1,ymm7                     ;distance+=(..)^2
    ;printregyps ymm1
    add r12,16                           ;avanzo di indicee

    jmp ciclo
resto:
    mov r11,[rbp+subb]                 ;subb
    sub r11,8                          ;subb-8
    ;vprintreg r11
ciclo2:
    cmp r12, r11
    jg fine
    vmovaps ymm0,[rdx+4*r12]             ;c[j*k*sub+i*sub+z]
    ;printregyps ymm0
    vsubps ymm0,[r8+4*r12]              ;c[j*k*sub+i*sub+z] - c[j*k*sub+j_d*sub+z]
    ;printregyps ymm0
    vmulps ymm0,ymm0                     ;(..)^2
    ;printregyps ymm0
    vaddps ymm2,ymm0                     ;distance+=(..)^2
    ;printregyps ymm2
    add r12,8                           ;avanzo di indice
    jmp ciclo2
fine:
    vaddps ymm1,ymm2        ;merge di tutte le somme
    ;printregyps ymm1
    vhaddps ymm1,ymm1        ;|
    ;printregyps ymm1
    vhaddps ymm1,ymm1        ;|
    ;printregyps ymm1
    vperm2f128 ymm5,ymm1,ymm1,1
    ;printregyps ymm5
    vaddss xmm1,xmm5
    ;printregyps ymm1
    vmovss  [rsi],xmm1
    stop

global rowDistance64SdcU

    rowDistance64SdcU:

    push		rbp				; salva il Base Pointer
    mov		rbp, rsp			; il Base Pointer punta al Record di Attivazione corrente
    pushaq						; salva i registri generali

    mov r11,[rbp+subb]

    
    ;mov dword [distance],0                    ;distance=0
      
    imul rcx,dim                        ;prendo j
    imul rcx,r11                        ;j*sub
    imul rcx,r9                        ;j*k*subb

    imul rdx,dim                        ;prendo i
    imul rdx,r11                        ;i*subb
    add rdx,rcx                         ;j*k*subb+i*subb
    add rdx,rdi                         ;c[j*k*subb+i*subb]
 
    imul r8,dim                         ;prendo j_d
    imul r8,r11                        ;j_d*subb
    add r8,rcx                         ;j*k*subb+j_d*subb
    add r8,rdi                         ;c[j*k*subb+j_d*subb]
    
    ;printreg rdx
    ;printreg r9
    vmovups ymm1,[rdx]                   ;c[j*k*sub+i*sub+0]
    vsubps ymm1,[r8]                    ;c[j*k*sub+i*sub+0] - c[j*k*sub+j_d*sub+0]
    vmulps ymm1,ymm1                     ;(..)^2
    ;printregyps ymm1
    ;mov rdi,dim
    sub r11,16                          ;subb-16

    mov r12,8                           ;z=4
    vxorps ymm2, ymm2
cicloU:
    cmp r12,r11
    jg restoU

    vmovups ymm0,[rdx+4*r12]             ;c[j*k*sub+i*sub+z]
    ;printregyps ymm0
    vmovups ymm7,[rdx+4*r12+32]             ;c[j*k*sub+i*sub+z]
    ;printregyps ymm0
    vsubps ymm0,[r8+4*r12]              ;c[j*k*sub+i*sub+z] - c[j*k*sub+j_d*sub+z]
    vsubps ymm7,[r8+4*r12+32]              ;c[j*k*sub+i*sub+z] - c[j*k*sub+j_d*sub+z]
    ;printregyps ymm0
    vmulps ymm0,ymm0                     ;(..)^2
    vmulps ymm7,ymm7                     ;(..)^2
    ;printregyps ymm0
    vaddps ymm1,ymm0                     ;distance+=(..)^2
    vaddps ymm1,ymm7                     ;distance+=(..)^2
    ;printregyps ymm1
    add r12,16                           ;avanzo di indice
    jmp cicloU
restoU:
    mov r11,[rbp+subb]                 ;subb
    sub r11,8                          ;subb-
    ;vprintreg r11
ciclo2U:
    cmp r12, r11
    jg resto2U
    vmovups ymm0,[rdx+4*r12]             ;c[j*k*sub+i*sub+z]
    ;printregyps ymm0
    vsubps ymm0,[r8+4*r12]              ;c[j*k*sub+i*sub+z] - c[j*k*sub+j_d*sub+z]
    ;printregyps ymm0
    vmulps ymm0,ymm0                     ;(..)^2
    ;printregyps ymm0
    vaddps ymm2,ymm0                     ;distance+=(..)^2
    ;printregyps ymm1
    add r12,8                           ;avanzo di indice
    jmp ciclo2U
resto2U:
    mov r11, [rbp+subb]
ciclo3U:
    cmp r12, r11
    je  fineU
    vmovss xmm7,[rdx+4*r12] ;sommo gli ultimi elementi rimanenti
    vsubss xmm7,[r8+4*r12]
    vmulss xmm7,xmm7
    vaddps ymm2, ymm7
    ;printregyps ymm2
    add r12, 1
    jmp ciclo3U
fineU:    
    vaddps ymm1,ymm2        ;merge di tutte le somme
    ;printregyps ymm1
    vhaddps ymm1,ymm1        ;|
    ;printregyps ymm1
    vhaddps ymm1,ymm1        ;|
    ;printregyps ymm1
    vperm2f128 ymm5,ymm1,ymm1,1
    vaddss xmm1,xmm5
    ;printregyps ymm1
    vmovss  [rsi],xmm1
    stop
