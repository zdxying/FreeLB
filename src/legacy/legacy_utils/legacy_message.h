// message manager
// legacy_message.h
#pragma once

#include <iostream>
#include <string>
#include <vector>

struct message {
  std::vector<bool> msg;

  message(int num) : msg(num) {
    for (int i = 0; i < num; i++) {
      msg[i] = false;
    }
  }

  void msg_out(int i, std::string out) {
    if (!msg[i]) {
      std::cout << out << std::endl;
      msg[i] = true;
    }
  }
};
