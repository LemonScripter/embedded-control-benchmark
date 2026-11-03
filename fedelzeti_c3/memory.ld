ENTRY(_start)
MEMORY {
  FLASH (rx)  : ORIGIN = 0x42000000, LENGTH = 4M
  RAM   (rwx) : ORIGIN = 0x3FC80000, LENGTH = 384K
}
SECTIONS {
  .text : { KEEP(*(.text.init)) *(.text .text.*) } > FLASH
  .rodata : { *(.rodata .rodata.*) } > FLASH
  .data : { *(.data .data.*) } > RAM
  .bss (NOLOAD) : { *(.bss .bss.*) *(COMMON) } > RAM
  _stack_top = ORIGIN(RAM) + LENGTH(RAM);
  /DISCARD/ : { *(.eh_frame*) *(.comment) }
}
